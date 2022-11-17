# Code for coordinating events on the printer toolhead
#
# Copyright (C) 2016-2021  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import math, logging, importlib, traceback, time
import mcu, chelper, kinematics.extruder
MoveType=chelper.MoveType
# Common suffixes: _d is distance (in mm), _v is velocity (in
#   mm/second), _v2 is velocity squared (mm^2/s^2), _t is time (in
#   seconds), _r is ratio (scalar between 0.0 and 1.0)

# Class to track each move request
# Note when jerk is not None, this is not a trapezoidal move. It is just a single accel segment.
class Move:
    def __init__(self, toolhead, start_pos, end_pos, speed,
                 secs, start_accel, jerk,
                 ext_end_v, ext_start_accel, ext_jerk):
        self.toolhead = toolhead
        self.start_pos = tuple(start_pos)
        self.end_pos = tuple(end_pos)
        self.axes_d = axes_d = [end_pos[i] - start_pos[i] for i in (0, 1, 2, 3)]
        self.move_d = move_d = math.sqrt(sum([d*d for d in axes_d[:3]]))
        if abs(self.axes_d[0] - -0.37713) < 1e-4:
            assert jerk is not None
        if jerk is not None:
            self.moveType = MoveType.withJerk
        else:
            self.moveType = MoveType.trapezoidal
        # accel aka start_accel
        self.secs = secs
        self.end_v = speed
        self.start_accel = start_accel
        self.accel = start_accel if start_accel is not None else toolhead.max_accel
        self.jerk = jerk
        self.ext_start_accel = ext_start_accel
        self.ext_jerk = ext_jerk
        self.ext_end_v = ext_end_v
        if self.moveType == MoveType.withJerk:
            print("fuck G7", self.start_pos, self.end_pos, self.axes_d, secs, speed, start_accel, jerk)
            #end_accel = start_accel + jerk * secs
            self.start_v = speed - start_accel * secs - 0.5*jerk*(secs**2)
            #self.start_v = speed - secs*(start_accel + 0.5*jerk*secs)
            self.ext_start_v = ext_end_v - ext_start_accel * secs - 0.5*ext_jerk*(secs**2)
            #assert self.start_v >= 0 and self.ext_start_v >= 0 and self.end_v >= 0
            # self.ext_start_v can be negative.
            assert self.start_v >= -1e-5 and self.end_v >= -1e-5, "%g %g" % (self.start_v, self.end_v)
            self.start_v = max(0, self.start_v)
            self.end_v = max(0, self.end_v)
            if self.start_v == 0:
                assert abs(self.ext_start_v) < 1e-5, "%g" % (self.ext_start_v)
                self.ext_start_v = 0
            if self.end_v == 0:
                assert self.ext_end_v == 0, "%g" % (self.ext_end_v)
            #print('fuck', self.start_v, secs, start_accel, jerk)
            expected_move_d = self.start_v * secs + 0.5 * start_accel * (secs**2) + 1.0/6.0*jerk*(secs**3)
            # TODO - this check needs to depend on precision of source gcode
            if abs(expected_move_d - self.move_d) > 1e-3:
                raise self.move_error("move paramters not self consistent distance %.5g vs %.5g (%.5g)" % (
                    expected_move_d, self.move_d, expected_move_d - self.move_d))

            ext_dist = self.end_pos[3] - self.start_pos[3]
            expect_ext_dist = self.ext_start_v * secs + 0.5 * self.ext_start_accel * (secs**2) + 1.0/6.0*self.ext_jerk*(secs**3)
            if abs(ext_dist - expect_ext_dist) > 1e-4:
                raise self.move_error("move extruder paramters not self consistent distance %.5g vs %.5g (%.5g)" % (
                    expect_ext_dist, ext_dist, expect_ext_dist - ext_dist))
            
        else:
            assert secs == None and start_accel == None
        
        self.timing_callbacks = []
        velocity = min(speed, toolhead.max_velocity)
        self.is_kinematic_move = True
        if move_d < .000000001:
            # Extrude only move
            self.end_pos = (start_pos[0], start_pos[1], start_pos[2],
                            end_pos[3])
            axes_d[0] = axes_d[1] = axes_d[2] = 0.
            self.move_d = move_d = abs(axes_d[3])
            inv_move_d = 0.
            if move_d:
                inv_move_d = 1. / move_d
            self.accel = 99999999.9
            velocity = speed
            self.is_kinematic_move = False
        else:
            inv_move_d = 1. / move_d
        self.axes_r = [d * inv_move_d for d in axes_d]
        # Junction speeds are tracked in velocity squared.  The
        # delta_v2 is the maximum amount of this squared-velocity that
        # can change in this move.
        if self.moveType == MoveType.trapezoidal:
            self.max_start_v2 = 0.
            self.max_cruise_v2 = velocity**2
            self.max_smoothed_v2 = 0.
            self.delta_v2 = 2.0 * move_d * self.accel
            self.smooth_delta_v2 = 2.0 * move_d * toolhead.max_accel_to_decel
            self.min_move_t = move_d / velocity
        else:
            if speed > toolhead.max_velocity or speed > toolhead.max_accel_to_decel:
                raise self.move_error("Invalid G7, velocity limited by toolhead.max_velocity or max_accel_to_decel %.5g" % (
                    speed))
            self.max_start_v2 = self.max_cruise_v2 = self.max_smoothed_v2 = self.delta_v2 = self.smooth_delta_v2 = None
            self.min_move_t = secs
    def limit_speed(self, speed, accel):
        speed2 = speed**2
        if self.moveType == MoveType.withJerk:
            # Not a perfect check because jerk could mean that extremum is in the middle somewhere.
            # Note accel can be negative...
            if speed < self.start_v or speed < self.end_v or accel != self.accel:
                raise self.move_error("move with jerk violates some limit %.5g vs %.5g, %g and %.5g,%.5g" % (speed, self.start_v, self.end_v, accel, self.accel))
            return
        if speed2 < self.max_cruise_v2:
            self.max_cruise_v2 = speed2
            self.min_move_t = self.move_d / speed
        self.accel = min(self.accel, accel)
        self.delta_v2 = 2.0 * self.move_d * self.accel
        self.smooth_delta_v2 = min(self.smooth_delta_v2, self.delta_v2)
    def move_error(self, msg="Move out of range"):
        print(self.axes_d)
        traceback.print_stack()
        ep = self.end_pos
        m = "%s: %.3f %.3f %.3f [%.3f]" % (msg, ep[0], ep[1], ep[2], ep[3])
        return self.toolhead.printer.command_error(m)
    def calc_junction(self, prev_move):
        assert self.moveType == MoveType.trapezoidal
        if not self.is_kinematic_move or not prev_move.is_kinematic_move:
            return
        # Allow extruder to calculate its maximum junction
        extruder_v2 = self.toolhead.extruder.calc_junction(prev_move, self)
        # Find max velocity using "approximated centripetal velocity"
        axes_r = self.axes_r
        prev_axes_r = prev_move.axes_r
        junction_cos_theta = -(axes_r[0] * prev_axes_r[0]
                               + axes_r[1] * prev_axes_r[1]
                               + axes_r[2] * prev_axes_r[2])
        if junction_cos_theta > 0.999999:
            return
        junction_cos_theta = max(junction_cos_theta, -0.999999)
        sin_theta_d2 = math.sqrt(0.5*(1.0-junction_cos_theta))
        R = (self.toolhead.junction_deviation * sin_theta_d2
             / (1. - sin_theta_d2))
        # Approximated circle must contact moves no further away than mid-move
        tan_theta_d2 = sin_theta_d2 / math.sqrt(0.5*(1.0+junction_cos_theta))
        move_centripetal_v2 = .5 * self.move_d * tan_theta_d2 * self.accel
        prev_move_centripetal_v2 = (.5 * prev_move.move_d * tan_theta_d2
                                    * prev_move.accel)
        # Apply limits
        if prev_move.moveType == MoveType.withJerk:
            # TODO - constrain min v2 as well.
            self.max_start_v2 = self.max_smoothed_v2 = prev_move.end_v ** 2
            return
        self.max_start_v2 = min(
            R * self.accel, R * prev_move.accel,
            move_centripetal_v2, prev_move_centripetal_v2,
            extruder_v2, self.max_cruise_v2, prev_move.max_cruise_v2,
            prev_move.max_start_v2 + prev_move.delta_v2)
        self.max_smoothed_v2 = min(
            self.max_start_v2
            , prev_move.max_smoothed_v2 + prev_move.smooth_delta_v2)
    def set_junction(self, start_v2, cruise_v2, end_v2):
        assert self.moveType == MoveType.trapezoidal
        # Determine accel, cruise, and decel portions of the move distance
        half_inv_accel = .5 / self.accel
        accel_d = (cruise_v2 - start_v2) * half_inv_accel
        decel_d = (cruise_v2 - end_v2) * half_inv_accel
        cruise_d = self.move_d - accel_d - decel_d
        # Determine move velocities
        self.start_v = start_v = math.sqrt(start_v2)
        self.cruise_v = cruise_v = math.sqrt(cruise_v2)
        self.end_v = end_v = math.sqrt(end_v2)
        # Determine time spent in each portion of move (time is the
        # distance divided by average velocity)
        self.accel_t = accel_d / ((start_v + cruise_v) * 0.5)
        self.cruise_t = cruise_d / cruise_v
        self.decel_t = decel_d / ((end_v + cruise_v) * 0.5)
        ext_axis_r = self.axes_r[3]
        self.ext_start_v = start_v * ext_axis_r
        self.ext_end_v = end_v * ext_axis_r

LOOKAHEAD_FLUSH_TIME = 0.250

# Class to track a list of pending move requests and to facilitate
# "look-ahead" across moves to reduce acceleration between moves.
class MoveQueue:
    def __init__(self, toolhead):
        self.toolhead = toolhead
        self.queue = []
        self.junction_flush = LOOKAHEAD_FLUSH_TIME
        #self.last_end_v = 0
        #self.ext_last_end_v = 0
    def reset(self):
        del self.queue[:]
        self.junction_flush = LOOKAHEAD_FLUSH_TIME
    def set_flush_time(self, flush_time):
        self.junction_flush = flush_time
    def get_last(self):
        if self.queue:
            return self.queue[-1]
        return None
    
    def flush(self, lazy=False):
        self.junction_flush = LOOKAHEAD_FLUSH_TIME
        update_flush_count = lazy
        queue = self.queue
        flush_count = len(queue)
        # Traverse queue from last to first move and determine maximum
        # junction speed assuming the robot comes to a complete stop
        # after the last move.
        #
        # If update_flush_count, then skip over some number of moves
        # at the end of the queue before start processing.
        delayed = []
        next_end_v2 = next_smoothed_v2 = peak_cruise_v2 = 0.
        next_move = None
        next_start_v2 = 0.
        print("START flush of len=", len(queue), lazy)
        for i in range(flush_count-1, -1, -1):
            move = queue[i]
            print("  %d(%s): %s" % (
                i, "G7" if move.moveType == MoveType.withJerk else "G0",
                str(move.axes_d)))
            if move.moveType == MoveType.withJerk:
                print("    start_v=%g end_v=%g ext_start_v=%g ext_end_v=%g" % (
                    move.start_v, move.end_v, move.ext_start_v, move.ext_end_v))
                if delayed:
                    # next_start_v2 hasn't been updated to next move yet.
                    raise move.move_error("Delayed moves after G7 move")
                start_v2 = move.start_v**2
                next_end_v2 = start_v2
                next_start_v2 = start_v2
                next_smoothed_v2 = start_v2
                move.accel_t = move.secs
                move.cruise_t = 0
                move.decel_t = 0
                next_move = move
                if not update_flush_count and i == (len(self.queue) - 1):
                    assert move.end_v == 0
                continue
            reachable_start_v2 = next_end_v2 + move.delta_v2
            start_v2 = min(move.max_start_v2, reachable_start_v2)
            reachable_smoothed_v2 = next_smoothed_v2 + move.smooth_delta_v2
            smoothed_v2 = min(move.max_smoothed_v2, reachable_smoothed_v2)
            actual_start_v2 = 0.
            did_set = False
            print("    reach=%.5g next=%.5g delta=%.5g v=%.5g peakc=%.5g" % (
                  math.sqrt(reachable_smoothed_v2), math.sqrt(next_smoothed_v2),
                  math.sqrt(move.smooth_delta_v2), math.sqrt(smoothed_v2),
                math.sqrt(peak_cruise_v2)))
            if smoothed_v2 < reachable_smoothed_v2:
                # It's possible for this move to accelerate
                if (smoothed_v2 + move.smooth_delta_v2 > next_smoothed_v2
                    or delayed):
                    # This move can decelerate or this is a full accel
                    # move after a full decel move
                    if update_flush_count and peak_cruise_v2:
                        print("    Setting flush_count=", i)
                        flush_count = i
                        update_flush_count = False
                    peak_cruise_v2 = min(move.max_cruise_v2, (
                        smoothed_v2 + reachable_smoothed_v2) * .5)
                    if delayed:
                        # Propagate peak_cruise_v2 to any delayed moves
                        if not update_flush_count and i < flush_count:
                            mc_v2 = peak_cruise_v2
                            for m, ms_v2, me_v2 in reversed(delayed):
                                mc_v2 = min(mc_v2, ms_v2)
                                actual_start_v2 = min(ms_v2, mc_v2)
                                actual_end_v2   = min(me_v2, mc_v2)
                                m.set_junction(actual_start_v2, mc_v2
                                               , actual_end_v2)
                                print("    (delayed) set start_v=%.5g end_v=%.5g" % (
                                    math.sqrt(actual_start_v2), math.sqrt(actual_end_v2)))
                                next_move = m
                        del delayed[:]
                if not update_flush_count and i < flush_count:
                    cruise_v2 = min((start_v2 + reachable_start_v2) * .5
                                    , move.max_cruise_v2, peak_cruise_v2)
                    actual_start_v2 = min(start_v2, cruise_v2)
                    actual_end_v2   = min(next_end_v2, cruise_v2)
                    #if abs(actual_end_v2 - next_start_v2) > 1e-8:
                    #    raise move.move_error("Discontinuous veloctiy change between moves %.5g vs %.5g" % (
                    #        math.sqrt(actual_end_v2), math.sqrt(next_start_v2)))
                    did_set = True
                    move.set_junction(actual_start_v2, cruise_v2, actual_end_v2)
                    assert not delayed
                    if i == (len(self.queue) - 1):
                        assert move.end_v == 0
                    if len(self.queue) == 1:
                        assert actual_start_v2 == 0
                    print("    set start_v=%.5g end_v=%.5g" % (
                        math.sqrt(actual_start_v2), math.sqrt(actual_end_v2)))
                    next_move = move
            else:
                # Delay calculating this move until peak_cruise_v2 is known
                did_set = True
                print("     delay append", next_move.axes_d if next_move else "None")
                delayed.append((move, start_v2, next_end_v2))
            next_start_v2 = actual_start_v2
            next_end_v2 = start_v2
            next_smoothed_v2 = smoothed_v2
            next_move = move
            if not did_set:
                print("    Skipping")
        if update_flush_count or not flush_count:
            return
        # Generate step times for all moves ready to be flushed
        self.toolhead._process_moves(queue[:flush_count])
        # Remove processed moves from the queue
        del queue[:flush_count]
    def add_move(self, move):
        self.queue.append(move)
        if len(self.queue) == 1:
            return
        if move.moveType == MoveType.trapezoidal:
            move.calc_junction(self.queue[-2])
        self.junction_flush -= move.min_move_t
        print("add_move", self.junction_flush, move.min_move_t)
        if self.junction_flush <= 0.:
            # Enough moves have been queued to reach the target flush time.
            self.flush(lazy=True)

MIN_KIN_TIME = 0.100
MOVE_BATCH_TIME = 0.500
SDS_CHECK_TIME = 0.001 # step+dir+step filter in stepcompress.c

DRIP_SEGMENT_TIME = 0.050
DRIP_TIME = 0.100
class DripModeEndSignal(Exception):
    pass

# Main code to track events (and their timing) on the printer toolhead
class ToolHead:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.reactor = self.printer.get_reactor()
        self.all_mcus = [
            m for n, m in self.printer.lookup_objects(module='mcu')]
        self.mcu = self.all_mcus[0]
        self.can_pause = True
        if self.mcu.is_fileoutput():
            self.can_pause = False
        self.move_queue = MoveQueue(self)
        self.commanded_pos = [0., 0., 0., 0.]
        self.printer.register_event_handler("klippy:shutdown",
                                            self._handle_shutdown)
        # Velocity and acceleration control
        self.max_velocity = config.getfloat('max_velocity', above=0.)
        self.max_accel = config.getfloat('max_accel', above=0.)
        self.requested_accel_to_decel = config.getfloat(
            'max_accel_to_decel', self.max_accel * 0.5, above=0.)
        self.max_accel_to_decel = self.requested_accel_to_decel
        self.square_corner_velocity = config.getfloat(
            'square_corner_velocity', 5., minval=0.)
        self.junction_deviation = 0.
        self._calc_junction_deviation()
        # Print time tracking
        self.buffer_time_low = config.getfloat(
            'buffer_time_low', 1.000, above=0.)
        self.buffer_time_high = config.getfloat(
            'buffer_time_high', 2.000, above=self.buffer_time_low)
        self.buffer_time_start = config.getfloat(
            'buffer_time_start', 0.250, above=0.)
        self.move_flush_time = config.getfloat(
            'move_flush_time', 0.050, above=0.)
        self.print_time = 0.
        self.special_queuing_state = "Flushed"
        self.need_check_stall = -1.
        self.flush_timer = self.reactor.register_timer(self._flush_handler)
        self.move_queue.set_flush_time(self.buffer_time_high)
        self.idle_flush_print_time = 0.
        self.print_stall = 0
        self.drip_completion = None
        # Kinematic step generation scan window time tracking
        self.kin_flush_delay = SDS_CHECK_TIME
        self.kin_flush_times = []
        self.last_kin_flush_time = self.last_kin_move_time = 0.
        # Setup iterative solver
        ffi_main, ffi_lib = chelper.get_ffi()
        self.trapq = ffi_main.gc(ffi_lib.trapq_alloc(), ffi_lib.trapq_free)
        self.trapq_append = ffi_lib.trapq_append
        self.trapq_append2 = ffi_lib.trapq_append2
        self.trapq_finalize_moves = ffi_lib.trapq_finalize_moves
        self.step_generators = []
        self.last_end_v = [0, 0, None]
        # Create kinematics class
        gcode = self.printer.lookup_object('gcode')
        self.Coord = gcode.Coord
        self.extruder = kinematics.extruder.DummyExtruder(self.printer)
        kin_name = config.get('kinematics')
        try:
            mod = importlib.import_module('kinematics.' + kin_name)
            self.kin = mod.load_kinematics(self, config)
        except config.error as e:
            raise
        except self.printer.lookup_object('pins').error as e:
            raise
        except:
            msg = "Error loading kinematics '%s'" % (kin_name,)
            logging.exception(msg)
            raise config.error(msg)
        # Register commands
        gcode.register_command('G4', self.cmd_G4)
        gcode.register_command('M400', self.cmd_M400)
        gcode.register_command('SET_VELOCITY_LIMIT',
                               self.cmd_SET_VELOCITY_LIMIT,
                               desc=self.cmd_SET_VELOCITY_LIMIT_help)
        gcode.register_command('M204', self.cmd_M204)
        # Load some default modules
        modules = ["gcode_move", "homing", "idle_timeout", "statistics",
                   "manual_probe", "tuning_tower"]
        for module_name in modules:
            self.printer.load_object(config, module_name)
    # Print time tracking
    def _update_move_time(self, next_print_time):
        batch_time = MOVE_BATCH_TIME
        kin_flush_delay = self.kin_flush_delay
        lkft = self.last_kin_flush_time
        while 1:
            self.print_time = min(self.print_time + batch_time, next_print_time)
            sg_flush_time = max(lkft, self.print_time - kin_flush_delay)
            for sg in self.step_generators:
                sg(sg_flush_time)
            free_time = max(lkft, sg_flush_time - kin_flush_delay)
            self.trapq_finalize_moves(self.trapq, free_time)
            self.extruder.update_move_time(free_time)
            mcu_flush_time = max(lkft, sg_flush_time - self.move_flush_time)
            for m in self.all_mcus:
                m.flush_moves(mcu_flush_time)
            if self.print_time >= next_print_time:
                break
    def _calc_print_time(self):
        curtime = self.reactor.monotonic()
        est_print_time = self.mcu.estimated_print_time(curtime)
        kin_time = max(est_print_time + MIN_KIN_TIME, self.last_kin_flush_time)
        kin_time += self.kin_flush_delay
        min_print_time = max(est_print_time + self.buffer_time_start, kin_time)
        if min_print_time > self.print_time:
            self.print_time = min_print_time
            self.printer.send_event("toolhead:sync_print_time",
                                    curtime, est_print_time, self.print_time)

    def check_junction(self, move):
        #print("check_junction", move.start_pos, move.end_pos)
        if abs(self.last_end_v[0] - move.start_v) > 1e-4:
            raise move.move_error("end_v doesn't match start_v of next move %.8f vs %.8f (%.5g)" % (
                self.last_end_v[0], move.start_v, self.last_end_v[0] - move.start_v))
        if move.moveType == MoveType.withJerk:
            ext_start_v = move.ext_start_v
            ext_end_v = move.ext_end_v
        else:
            ext_start_v = move.start_v * move.axes_r[3]
            ext_end_v = move.end_v * move.axes_r[3]
        #print("check_junction: ", ext_start_v, ext_end_v, self.last_end_v)
        #print(move.move_d, move.axes_d[3] / move.move_d, move.axes_d[3] / move.move_d * move.start_v, move.start_v, move.axes_r[3], ext_start_v)
        if move.moveType == MoveType.withJerk and abs(self.last_end_v[1] - ext_start_v) > 1e-4:
            # If mvoe is trapezoidal, it could have an instantaneous extruder change
            # due to extruder instant_corner_v
            raise move.move_error("extruder end_v doesn't match start_v of next move %.5g vs %.5g" % (
                self.last_end_v[1], ext_start_v))

        tmove = move
        if self.special_queuing_state == "Drip":
            tmove = None
        elif self.last_end_v[2] is not None:
            for i in (0,1,2,3):
                if abs(self.last_end_v[2].end_pos[i] - move.start_pos[i]) > 1e-3:
                    print(self.special_queuing_state)
                    raise move.move_error("Illegal instantaneous position change in move %d: %g->%g" % (i, self.last_end_v[2].end_pos[i], move.start_pos[i]))
            
        self.last_end_v = [ move.end_v, ext_end_v, tmove ]

    def _process_moves(self, moves):
        # Resync print_time if necessary
        if self.special_queuing_state:
            if self.special_queuing_state != "Drip":
                # Transition from "Flushed"/"Priming" state to main state
                self.special_queuing_state = ""
                self.need_check_stall = -1.
                self.reactor.update_timer(self.flush_timer, self.reactor.NOW)
            self._calc_print_time()
        # Queue moves into trapezoid motion queue (trapq)
        next_move_time = self.print_time
        for move in moves:
            print("trapq_append(%s): %s %s  start=%.5g %send=%.5g" % (
                "G7" if move.moveType == MoveType.withJerk else "G0",
                str(move.axes_d), str(move.end_pos), move.start_v,
                ("cruise=%.5g "%move.cruise_v) if move.moveType == MoveType.trapezoidal else "",
                move.end_v))
            if move.is_kinematic_move:
                if move.moveType == MoveType.withJerk:
                    self.trapq_append2(
                        self.trapq, next_move_time,
                        move.secs,
                        move.start_pos[0], move.start_pos[1], move.start_pos[2],
                        move.axes_r[0], move.axes_r[1], move.axes_r[2],
                        move.start_v, move.start_accel, move.jerk)
                else:
                    self.trapq_append(
                        self.trapq, next_move_time,
                        move.accel_t, move.cruise_t, move.decel_t,
                        move.start_pos[0], move.start_pos[1], move.start_pos[2],
                        move.axes_r[0], move.axes_r[1], move.axes_r[2],
                        move.start_v, move.cruise_v, move.accel)
            self.check_junction(move)
            if move.axes_d[3]:
                self.extruder.move(next_move_time, move)
            next_move_time = (next_move_time + move.accel_t
                              + move.cruise_t + move.decel_t)
            for cb in move.timing_callbacks:
                cb(next_move_time)
        # Generate steps for moves
        if self.special_queuing_state:
            self._update_drip_move_time(next_move_time)
        self._update_move_time(next_move_time)
        self.last_kin_move_time = next_move_time
    def flush_step_generation(self):
        # Transition from "Flushed"/"Priming"/main state to "Flushed" state
        self.move_queue.flush()
        self.special_queuing_state = "Flushed"
        self.need_check_stall = -1.
        self.reactor.update_timer(self.flush_timer, self.reactor.NEVER)
        self.move_queue.set_flush_time(self.buffer_time_high)
        self.idle_flush_print_time = 0.
        flush_time = self.last_kin_move_time + self.kin_flush_delay
        flush_time = max(flush_time, self.print_time - self.kin_flush_delay)
        self.last_kin_flush_time = max(self.last_kin_flush_time, flush_time)
        self._update_move_time(max(self.print_time, self.last_kin_flush_time))
    def _flush_lookahead(self):
        if self.special_queuing_state:
            return self.flush_step_generation()
        self.move_queue.flush()
    def get_last_move_time(self):
        self._flush_lookahead()
        if self.special_queuing_state:
            self._calc_print_time()
        return self.print_time
    def _check_stall(self):
        eventtime = self.reactor.monotonic()
        if self.special_queuing_state:
            if self.idle_flush_print_time:
                # Was in "Flushed" state and got there from idle input
                est_print_time = self.mcu.estimated_print_time(eventtime)
                if est_print_time < self.idle_flush_print_time:
                    self.print_stall += 1
                self.idle_flush_print_time = 0.
            # Transition from "Flushed"/"Priming" state to "Priming" state
            self.special_queuing_state = "Priming"
            self.need_check_stall = -1.
            self.reactor.update_timer(self.flush_timer, eventtime + 0.100)
        # Check if there are lots of queued moves and stall if so
        while 1:
            est_print_time = self.mcu.estimated_print_time(eventtime)
            buffer_time = self.print_time - est_print_time
            stall_time = buffer_time - self.buffer_time_high
            if stall_time <= 0.:
                break
            if not self.can_pause:
                self.need_check_stall = self.reactor.NEVER
                return
            eventtime = self.reactor.pause(eventtime + min(1., stall_time))
        if not self.special_queuing_state:
            # In main state - defer stall checking until needed
            self.need_check_stall = (est_print_time + self.buffer_time_high
                                     + 0.100)
    def _flush_handler(self, eventtime):
        try:
            print_time = self.print_time
            buffer_time = print_time - self.mcu.estimated_print_time(eventtime)
            if buffer_time > self.buffer_time_low:
                # Running normally - reschedule check
                return eventtime + buffer_time - self.buffer_time_low
            # Under ran low buffer mark - flush lookahead queue
            self.flush_step_generation()
            if print_time != self.print_time:
                self.idle_flush_print_time = self.print_time
        except:
            logging.exception("Exception in flush_handler")
            self.printer.invoke_shutdown("Exception in flush_handler")
        return self.reactor.NEVER
    # Movement commands
    def get_position(self):
        return list(self.commanded_pos)
    def set_position(self, newpos, homing_axes=()):
        self.flush_step_generation()
        ffi_main, ffi_lib = chelper.get_ffi()
        ffi_lib.trapq_set_position(self.trapq, self.print_time,
                                   newpos[0], newpos[1], newpos[2])
        self.last_end_v[2] = None
        self.commanded_pos[:] = newpos
        self.kin.set_position(newpos, homing_axes)
        self.printer.send_event("toolhead:set_position")
    def move(self, newpos, speed,
             secs=None, start_accel=None, jerk=None,
             ext_end_v=None, ext_start_accel=None, ext_jerk=None):
        move = Move(self, self.commanded_pos, newpos, speed,
                    secs, start_accel, jerk,
                    ext_end_v, ext_start_accel, ext_jerk)
        #print("Move requested: ", newpos, speed, move.axes_d)
        #traceback.print_stack()
        if not move.move_d:
            return
        if move.is_kinematic_move:
            self.kin.check_move(move)
        if move.axes_d[3]:
            self.extruder.check_move(move)
        self.commanded_pos[:] = move.end_pos

        self.move_queue.add_move(move)
        if self.print_time > self.need_check_stall:
            self._check_stall()
    # Is print head halt-able at this point.  Invariant which G7 breaks
    # is that the MoveQueue always can end with v=0
    def last_move_can_pause(self):
        last_move = self.move_queue.get_last()
        return not last_move or \
            last_move.moveType == MoveType.trapezoidal or \
            last_move.end_v == 0
    def manual_move(self, coord, speed):
        curpos = list(self.commanded_pos)
        for i in range(len(coord)):
            if coord[i] is not None:
                curpos[i] = coord[i]
        self.move(curpos, speed)
        self.printer.send_event("toolhead:manual_move")
    def dwell(self, delay):
        print("Josh - sleeping for a second", delay)
        time.sleep(delay / 1000.)
        next_print_time = self.get_last_move_time() + max(0., delay)
        self._update_move_time(next_print_time)
        self._check_stall()
    def wait_moves(self):
        self._flush_lookahead()
        eventtime = self.reactor.monotonic()
        while (not self.special_queuing_state
               or self.print_time >= self.mcu.estimated_print_time(eventtime)):
            if not self.can_pause:
                break
            eventtime = self.reactor.pause(eventtime + 0.100)
    def set_extruder(self, extruder, extrude_pos):
        self.extruder = extruder
        self.commanded_pos[3] = extrude_pos
        self.last_end_v[2] = None
    def get_extruder(self):
        return self.extruder
    # Homing "drip move" handling
    def _update_drip_move_time(self, next_print_time):
        flush_delay = DRIP_TIME + self.move_flush_time + self.kin_flush_delay
        while self.print_time < next_print_time:
            if self.drip_completion.test():
                raise DripModeEndSignal()
            curtime = self.reactor.monotonic()
            est_print_time = self.mcu.estimated_print_time(curtime)
            wait_time = self.print_time - est_print_time - flush_delay
            if wait_time > 0. and self.can_pause:
                # Pause before sending more steps
                self.drip_completion.wait(curtime + wait_time)
                continue
            npt = min(self.print_time + DRIP_SEGMENT_TIME, next_print_time)
            self._update_move_time(npt)
    def drip_move(self, newpos, speed, drip_completion):
        self.dwell(self.kin_flush_delay)
        # Transition from "Flushed"/"Priming"/main state to "Drip" state
        self.move_queue.flush()
        self.special_queuing_state = "Drip"
        self.need_check_stall = self.reactor.NEVER
        self.reactor.update_timer(self.flush_timer, self.reactor.NEVER)
        self.move_queue.set_flush_time(self.buffer_time_high)
        self.idle_flush_print_time = 0.
        self.drip_completion = drip_completion
        # Submit move
        try:
            self.move(newpos, speed)
        except self.printer.command_error as e:
            self.flush_step_generation()
            raise
        # Transmit move in "drip" mode
        try:
            self.move_queue.flush()
        except DripModeEndSignal as e:
            self.move_queue.reset()
            self.trapq_finalize_moves(self.trapq, self.reactor.NEVER)
        # Exit "Drip" state
        self.flush_step_generation()
    # Misc commands
    def stats(self, eventtime):
        for m in self.all_mcus:
            m.check_active(self.print_time, eventtime)
        buffer_time = self.print_time - self.mcu.estimated_print_time(eventtime)
        is_active = buffer_time > -60. or not self.special_queuing_state
        if self.special_queuing_state == "Drip":
            buffer_time = 0.
        return is_active, "print_time=%.3f buffer_time=%.3f print_stall=%d" % (
            self.print_time, max(buffer_time, 0.), self.print_stall)
    def check_busy(self, eventtime):
        est_print_time = self.mcu.estimated_print_time(eventtime)
        lookahead_empty = not self.move_queue.queue
        return self.print_time, est_print_time, lookahead_empty
    def get_status(self, eventtime):
        print_time = self.print_time
        estimated_print_time = self.mcu.estimated_print_time(eventtime)
        res = dict(self.kin.get_status(eventtime))
        res.update({ 'print_time': print_time,
                     'stalls': self.print_stall,
                     'estimated_print_time': estimated_print_time,
                     'extruder': self.extruder.get_name(),
                     'position': self.Coord(*self.commanded_pos),
                     'max_velocity': self.max_velocity,
                     'max_accel': self.max_accel,
                     'max_accel_to_decel': self.requested_accel_to_decel,
                     'square_corner_velocity': self.square_corner_velocity})
        return res
    def _handle_shutdown(self):
        self.can_pause = False
        self.move_queue.reset()
    def get_kinematics(self):
        return self.kin
    def get_trapq(self):
        return self.trapq
    def register_step_generator(self, handler):
        self.step_generators.append(handler)
    def note_step_generation_scan_time(self, delay, old_delay=0.):
        self.flush_step_generation()
        cur_delay = self.kin_flush_delay
        if old_delay:
            self.kin_flush_times.pop(self.kin_flush_times.index(old_delay))
        if delay:
            self.kin_flush_times.append(delay)
        new_delay = max(self.kin_flush_times + [SDS_CHECK_TIME])
        self.kin_flush_delay = new_delay
    def register_lookahead_callback(self, callback):
        last_move = self.move_queue.get_last()
        if last_move is None:
            callback(self.get_last_move_time())
            return
        last_move.timing_callbacks.append(callback)
    def note_kinematic_activity(self, kin_time):
        self.last_kin_move_time = max(self.last_kin_move_time, kin_time)
    def get_max_velocity(self):
        return self.max_velocity, self.max_accel
    def _calc_junction_deviation(self):
        scv2 = self.square_corner_velocity**2
        self.junction_deviation = scv2 * (math.sqrt(2.) - 1.) / self.max_accel
        self.max_accel_to_decel = min(self.requested_accel_to_decel,
                                      self.max_accel)
    def cmd_G4(self, gcmd):
        # Dwell
        delay = gcmd.get_float('P', 0., minval=0.) / 1000.
        self.dwell(delay)
    def cmd_M400(self, gcmd):
        # Wait for current moves to finish
        self.wait_moves()
    cmd_SET_VELOCITY_LIMIT_help = "Set printer velocity limits"
    def cmd_SET_VELOCITY_LIMIT(self, gcmd):
        max_velocity = gcmd.get_float('VELOCITY', None, above=0.)
        max_accel = gcmd.get_float('ACCEL', None, above=0.)
        square_corner_velocity = gcmd.get_float(
            'SQUARE_CORNER_VELOCITY', None, minval=0.)
        requested_accel_to_decel = gcmd.get_float(
            'ACCEL_TO_DECEL', None, above=0.)
        if max_velocity is not None:
            self.max_velocity = max_velocity
        if max_accel is not None:
            self.max_accel = max_accel
        if square_corner_velocity is not None:
            self.square_corner_velocity = square_corner_velocity
        if requested_accel_to_decel is not None:
            self.requested_accel_to_decel = requested_accel_to_decel
        self._calc_junction_deviation()
        msg = ("max_velocity: %.6f\n"
               "max_accel: %.6f\n"
               "max_accel_to_decel: %.6f\n"
               "square_corner_velocity: %.6f" % (
                   self.max_velocity, self.max_accel,
                   self.requested_accel_to_decel,
                   self.square_corner_velocity))
        self.printer.set_rollover_info("toolhead", "toolhead: %s" % (msg,))
        if (max_velocity is None and
            max_accel is None and
            square_corner_velocity is None and
            requested_accel_to_decel is None):
            gcmd.respond_info(msg, log=False)
    def cmd_M204(self, gcmd):
        # Use S for accel
        accel = gcmd.get_float('S', None, above=0.)
        if accel is None:
            # Use minimum of P and T for accel
            p = gcmd.get_float('P', None, above=0.)
            t = gcmd.get_float('T', None, above=0.)
            if p is None or t is None:
                gcmd.respond_info('Invalid M204 command "%s"'
                                  % (gcmd.get_commandline(),))
                return
            accel = min(p, t)
        self.max_accel = accel
        self._calc_junction_deviation()

def add_printer_objects(config):
    config.get_printer().add_object('toolhead', ToolHead(config))
    kinematics.extruder.add_printer_objects(config)
