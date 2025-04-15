# Perform Z Homing at specific XY coordinates.
#
# Copyright (C) 2019 Florian Heilmann <Florian.Heilmann@gmx.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

class SafeZHoming:
    def __init__(self, config):
        self.printer = config.get_printer()
        x_pos, y_pos = config.getfloatlist("home_xy_position", count=2)
        self.home_x_pos, self.home_y_pos = x_pos, y_pos
        self.z_hop = config.getfloat("z_hop", default=0.0)
        self.z_hop_speed = config.getfloat('z_hop_speed', 15., above=0.)
        zconfig = config.getsection('stepper_z')
        self.max_z = zconfig.getfloat('position_max', note_valid=False)
        self.speed = config.getfloat('speed', 50.0, above=0.)
        self.move_to_previous = config.getboolean('move_to_previous', False)
        self.printer.load_object(config, 'homing')
        self.gcode = self.printer.lookup_object('gcode')
        self.prev_G28 = self.gcode.register_command("G28", None)
        self.gcode.register_command("G28", self.cmd_G28)
        self.gcode.register_command('JPROBE', self.cmd_JProbe)

        if config.has_section("homing_override"):
            raise config.error("homing_override and safe_z_homing cannot"
                               +" be used simultaneously")

    def cmd_JProbe(self, gcmd):
        toolhead = self.printer.lookup_object('toolhead')
        pprobe = self.printer.lookup_object('probe')
        x_offset, y_offset, z_offset = pprobe.get_offsets()

        # Home Z
        probe_session = pprobe.start_probe_session(gcmd)
        self.gcode.run_script_from_command('G28 Z')
        pos = toolhead.get_position()
        toolhead.manual_move([None, None, pos[2] + probe_session.sample_retract_dist], probe_session.speed)
        # Now do more detailed probing
        probe_session.run_probe(gcmd)
        rez = probe_session.pull_probed_results()
        assert len(rez) == 1, rez # Each run_probe() generates a single result
        assert len(rez[0]) == 3, rez[0]
        avgZpos = rez[0][2]
        
        pos = toolhead.get_position()
        # I think pos[2] reflects the last probe, but avgZpos reflects the average
        self.gcode.respond_info("jprobe: head z move since probe: %g (%g vs %g)" % (
            pos[2] - avgZpos, pos[2], avgZpos))
        #logging.debug("jprobe: ", rez, probe.z_offset, pos[2])
        if abs(pos[2] - avgZpos) > 1.0:
            raise gcmd.error("Somehow head moved after probing: %g vs %g" % (pos[2], avgZpos))
        pos[2] = z_offset + (pos[2] - avgZpos)
        assert pos[2] >= 0
        toolhead.set_position(pos, homing_axes=[2])
        toolhead.manual_move([None, None, pos[2] + probe_session.sample_retract_dist], self.z_hop_speed)
        probe_session.end_probe_session()
        
    def cmd_G28(self, gcmd):
        toolhead = self.printer.lookup_object('toolhead')

        # Perform Z Hop if necessary
        if self.z_hop != 0.0:
            # Check if Z axis is homed and its last known position
            curtime = self.printer.get_reactor().monotonic()
            kin_status = toolhead.get_kinematics().get_status(curtime)
            pos = toolhead.get_position()

            if 'z' not in kin_status['homed_axes']:
                # Always perform the z_hop if the Z axis is not homed
                pos[2] = 0
                toolhead.set_position(pos, homing_axes="z")
                toolhead.manual_move([None, None, self.z_hop],
                                     self.z_hop_speed)
                toolhead.get_kinematics().clear_homing_state("z")
            elif pos[2] < self.z_hop:
                # If the Z axis is homed, and below z_hop, lift it to z_hop
                toolhead.manual_move([None, None, self.z_hop],
                                     self.z_hop_speed)

        # Determine which axes we need to home
        need_x, need_y, need_z = [gcmd.get(axis, None) is not None
                                  for axis in "XYZ"]
        if not need_x and not need_y and not need_z:
            need_x = need_y = need_z = True

        # Home XY axes if necessary
        new_params = {}
        if need_x:
            new_params['X'] = '0'
        if need_y:
            new_params['Y'] = '0'
        if new_params:
            g28_gcmd = self.gcode.create_gcode_command("G28", "G28", new_params)
            self.prev_G28(g28_gcmd)

        # Home Z axis if necessary
        if need_z:
            # Throw an error if X or Y are not homed
            curtime = self.printer.get_reactor().monotonic()
            kin_status = toolhead.get_kinematics().get_status(curtime)
            if ('x' not in kin_status['homed_axes'] or
                'y' not in kin_status['homed_axes']):
                raise gcmd.error("Must home X and Y axes first")
            # Move to safe XY homing position
            prevpos = toolhead.get_position()
            toolhead.manual_move([self.home_x_pos, self.home_y_pos], self.speed)
            # Home Z
            g28_gcmd = self.gcode.create_gcode_command("G28", "G28", {'Z': '0'})
            self.prev_G28(g28_gcmd)
            # Perform Z Hop again for pressure-based probes
            if self.z_hop:
                pos = toolhead.get_position()
                if pos[2] < self.z_hop:
                    toolhead.manual_move([None, None, self.z_hop],
                                         self.z_hop_speed)
            # Move XY back to previous positions
            if self.move_to_previous:
                toolhead.manual_move(prevpos[:2], self.speed)

def load_config(config):
    return SafeZHoming(config)
