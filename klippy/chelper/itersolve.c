// Iterative solver for kinematic moves
//
// Copyright (C) 2018-2020  Kevin O'Connor <kevin@koconnor.net>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // fabs
#include <stddef.h> // offsetof
#include <string.h> // memset
#include "compiler.h" // __visible
#include "itersolve.h" // itersolve_generate_steps
#include "pyhelper.h" // errorf
#include "stepcompress.h" // queue_append_start
#include "trapq.h" // struct move
#include <stdio.h>

/****************************************************************
 * Main iterative solver
 ****************************************************************/

struct timepos {
    double time, position;
};

#define SEEK_TIME_RESET 0.000100

//#define DEUBG

#if 1
struct stepcompress {
    // Buffer management
    uint32_t *queue, *queue_end, *queue_pos, *queue_next;
    // Internal tracking
    uint32_t max_error;
    double mcu_time_offset, mcu_freq, last_step_print_time;
    // Message generation
    uint64_t last_step_clock;
    struct list_head msg_queue;
    uint32_t oid;
    int32_t queue_step_msgtag, set_next_step_dir_msgtag;
    int sdir, invert_sdir;
    // Step+dir+step filter
    uint64_t next_step_clock;
    int next_step_dir;
    // History tracking
    int64_t last_position;
    struct list_head history_list;
};
#endif

static inline int check_active(struct stepper_kinematics *sk, struct move *m);

// Generate step times for a portion of a move
static int32_t
itersolve_gen_steps_range(struct stepper_kinematics *sk, struct move *m
                          , double abs_start, double abs_end)
{
    sk_calc_callback calc_position_cb = sk->calc_position_cb;
    double half_step = .5 * sk->step_dist;
    double start = abs_start - m->print_time, end = abs_end - m->print_time;
    if (start < 0.)
        start = 0.;
    if (end > m->move_t)
        end = m->move_t;
    //printf("  gen steps2: start=%g end=%g\n", start, end);
    struct timepos old_guess = {start, sk->commanded_pos}, guess = old_guess;
    int sdir = stepcompress_get_step_dir(sk->sc);
    int is_dir_change = 0, have_bracket = 0, check_oscillate = 0;
    double target = sk->commanded_pos + (sdir ? half_step : -half_step);

    // could be invoking shaper_xy_calc_position...
    double t = calc_position_cb(sk, m, start);
#ifdef DEBUG
    printf("    do gen steps %d: t=%g[%g->%g] start_v=%g accel=%g jerk=%g start_pos=%g,%g,%g commanded_pos=%g astartpos=%g\n", sk->sc->oid,
           m->move_t, start, end, m->start_v, m->half_accel*2.0, m->sixth_jerk*6.0,
           m->start_pos.x, m->start_pos.y, m->start_pos.z, sk->commanded_pos, t);
#endif
    if (end <= 0) {
#ifdef DEBUG
      printf("       SKIPPING cuz negative end time\n");
#endif
      return 0;
    }
    if (end <= start) {
#ifdef DEBUG
      printf("       SKIPPING cuz end before start\n");
#endif
      return 0;
    }
    double instant_steps = fabs(t - sk->commanded_pos) / half_step;
    if (instant_steps > 4.0) {
      // We're asking for instantaneous jump - I don't see how that could work.
#ifdef DEBUG
      printf("      j-itersolve_gen_steps_range %d: Instantaneous position change by %g steps: move.start=%g commanded_pos=%g half_step=%g\n",
             sk->sc->oid, instant_steps, t, sk->commanded_pos, half_step);
#endif
      errorf("itersolve_gen_steps_range %d: Instantaneous position change by %g steps: move.start=%g commanded_pos=%g half_step=%g\n",
             sk->sc->oid, instant_steps, t, sk->commanded_pos, half_step);
      return ERROR_RET;
    }
    double last_time=start, low_time=start, high_time=start + SEEK_TIME_RESET;
    if (high_time > end)
        high_time = end;
    int stepCount = 0;
    for (;;) {
        // Use the "secant method" to guess a new time from previous guesses
        double guess_dist = guess.position - target;
        double og_dist = old_guess.position - target;
        double next_time = ((old_guess.time*guess_dist - guess.time*og_dist)
                            / (guess_dist - og_dist));
        if (!(next_time > low_time && next_time < high_time)) { // or NaN
            // Next guess is outside bounds checks - validate it
            if (have_bracket) {
                // A poor guess - fall back to bisection
                next_time = (low_time + high_time) * .5;
                check_oscillate = 0;
            } else if (guess.time >= end) {
                // No more steps present in requested time range
                break;
            } else {
                // Might be a poor guess - limit to exponential search
                next_time = high_time;
                high_time = 2. * high_time - last_time;
                if (high_time > end)
                    high_time = end;
            }
        }
        // Calculate position at next_time guess
        old_guess = guess;
        guess.time = next_time;
        guess.position = calc_position_cb(sk, m, next_time);
        guess_dist = guess.position - target;
        if (fabs(guess_dist) > .000000001) {
            // Guess does not look close enough - update bounds
            double rel_dist = sdir ? guess_dist : -guess_dist;
            if (rel_dist > 0.) {
                // Found position past target, so step is definitely present
                if (have_bracket && old_guess.time <= low_time) {
                    if (check_oscillate)
                        // Force bisect next to avoid persistent oscillations
                        old_guess = guess;
                    check_oscillate = 1;
                }
                high_time = guess.time;
                have_bracket = 1;
            } else if (rel_dist < -(half_step + half_step + .000000010)) {
                // Found direction change
                sdir = !sdir;
                target = (sdir ? target + half_step + half_step
                          : target - half_step - half_step);
                low_time = last_time;
                high_time = guess.time;
                is_dir_change = have_bracket = 1;
                check_oscillate = 0;
            } else {
                low_time = guess.time;
            }
            if (!have_bracket || high_time - low_time > .000000001) {
                if (!is_dir_change && rel_dist >= -half_step)
                    // Avoid rollback if stepper fully reaches step position
                    stepcompress_commit(sk->sc);
                // Guess is not close enough - guess again with new time
                continue;
            }
        }
        // Found next step - submit it
        //printf("  stepcompress_append: %g %g %g %g %g %g\n",
        //       m->print_time, guess.time, guess.position,
        //       guess_dist, target, half_step);
        stepCount++;
        int ret = stepcompress_append(sk->sc, sdir, m->print_time, guess.time);
        if (ret)
            return ret;
        target = sdir ? target+half_step+half_step : target-half_step-half_step;
        // Reset bounds checking
        double seek_time_delta = 1.5 * (guess.time - last_time);
        if (seek_time_delta < .000000001)
            seek_time_delta = .000000001;
        if (is_dir_change && seek_time_delta > SEEK_TIME_RESET)
            seek_time_delta = SEEK_TIME_RESET;
        last_time = low_time = guess.time;
        high_time = guess.time + seek_time_delta;
        if (high_time > end)
            high_time = end;
        is_dir_change = have_bracket = check_oscillate = 0;
    }
    sk->commanded_pos = target - (sdir ? half_step : -half_step);
    if (!check_active(sk, m) && stepCount > 0) {
#ifdef DEBUG
      printf("      WARNING j-itersolve_gen_steps_range %d: inactive move generated stepCount=%d start_pos=%g,%g,%g pt=%g mt=%g %g %g %g\n",
             sk->sc->oid, stepCount,
             m->start_pos.x, m->start_pos.y, m->start_pos.z,
             m->print_time, m->move_t, m->start_v, m->half_accel, m->sixth_jerk
             );
#endif
    }
    if (sk->post_cb)
        sk->post_cb(sk);
    return 0;
}


/****************************************************************
 * Interface functions
 ****************************************************************/

// Check if a move is likely to cause movement on a stepper
static inline int
check_active(struct stepper_kinematics *sk, struct move *m)
{
    int af = sk->active_flags;
    return ((af & AF_X && m->axes_r.x != 0.)
            || (af & AF_Y && m->axes_r.y != 0.)
            || (af & AF_Z && m->axes_r.z != 0.));
}

// Generate step times for a range of moves on the trapq
int32_t __visible
itersolve_generate_steps(struct stepper_kinematics *sk, double flush_time)
{
  if (flush_time < sk->last_flush_time) {
    printf("WARNING: oid=%d flush_time going backwards %g vs %g\n",
           sk->sc->oid,
           flush_time, sk->last_flush_time);
  }
#ifdef DEBUG
  printf("%d: start gen steps flush times %g -> %g\n", sk->sc->oid, sk->last_flush_time, flush_time);
#endif
  double last_flush_time = sk->last_flush_time;
    sk->last_flush_time = flush_time;
    if (!sk->tq)
        return 0;
    trapq_check_sentinels(sk->tq);
    struct move *m = list_first_entry(&sk->tq->moves, struct move, node);
    while (last_flush_time >= m->print_time + m->move_t)
        m = list_next_entry(m, node);
    double force_steps_time = sk->last_move_time + sk->gen_steps_post_active;
    int skip_count = 0;
    for (;;) {
#ifdef DEBUG
      printf("gen steps %d: pt=%g t=%g start_v=%g accel=%g jerk=%g start_pos=%g,%g,%g axes_r=%g,%g,%g\n", sk->sc->oid,
             m->print_time, m->move_t, m->start_v, m->half_accel*2.0, m->sixth_jerk*6.0,
             m->start_pos.x, m->start_pos.y, m->start_pos.z, m->axes_r.x, m->axes_r.y, m->axes_r.z);
#endif
      double move_start = m->print_time, move_end = move_start + m->move_t;
      //if (move_start >= flush_time) {
      //  printf("WARNING: oid=%d examining move in the future %g (%g vs %g)\n",
      //         sk->sc->oid,
      //         move_start - flush_time,
      //         flush_time, move_start);
      //}
        if (check_active(sk, m)) {
            if (skip_count && sk->gen_steps_pre_active) {
                // Must generate steps leading up to stepper activity
                double abs_start = move_start - sk->gen_steps_pre_active;
                if (abs_start < last_flush_time)
                    abs_start = last_flush_time;
                if (abs_start < force_steps_time)
                    abs_start = force_steps_time;
                struct move *pm = list_prev_entry(m, node);
                while (--skip_count && pm->print_time > abs_start)
                    pm = list_prev_entry(pm, node);
                do {
#ifdef DEBUG
                  printf("  gen steps prev move %d: pt=%g t=%g start_v=%g accel=%g jerk=%g start_pos=%g,%g\n", sk->sc->oid,
                         pm->print_time, pm->move_t, pm->start_v, pm->half_accel*2.0, pm->sixth_jerk*6.0,
                         pm->start_pos.x, pm->start_pos.y);
#endif
                  int32_t ret = itersolve_gen_steps_range(sk, pm, abs_start
                                                            , flush_time);
                    if (ret)
                        return ret;
                    pm = list_next_entry(pm, node);
                } while (pm != m);
            }
            // Generate steps for this move
#ifdef DEBUG
            printf("  gen steps this move %d t %g->%g\n", sk->sc->oid, last_flush_time, flush_time);
#endif
            int32_t ret = itersolve_gen_steps_range(sk, m, last_flush_time
                                                    , flush_time);
            if (ret)
                return ret;
            if (move_end >= flush_time) {
                sk->last_move_time = flush_time;
                return 0;
            }
            skip_count = 0;
            sk->last_move_time = move_end;
            force_steps_time = sk->last_move_time + sk->gen_steps_post_active;
        } else {
            if (move_start < force_steps_time) {
                // Must generates steps just past stepper activity
                double abs_end = force_steps_time;
                if (abs_end > flush_time)
                    abs_end = flush_time;
#ifdef DEBUG
                printf("  gen steps post move %d\n", sk->sc->oid);
#endif
                int32_t ret = itersolve_gen_steps_range(sk, m, last_flush_time
                                                        , abs_end);
                if (ret)
                    return ret;
                skip_count = 1;
            } else {
                // This move doesn't impact this stepper - skip it
                skip_count++;
            }
            if (flush_time + sk->gen_steps_pre_active <= move_end)
                return 0;
        }
        m = list_next_entry(m, node);
    }
}

// Check if the given stepper is likely to be active in the given time range
double __visible
itersolve_check_active(struct stepper_kinematics *sk, double flush_time)
{
    if (!sk->tq)
        return 0.;
    trapq_check_sentinels(sk->tq);
    struct move *m = list_first_entry(&sk->tq->moves, struct move, node);
    while (sk->last_flush_time >= m->print_time + m->move_t)
        m = list_next_entry(m, node);
    for (;;) {
        if (check_active(sk, m))
            return m->print_time;
        if (flush_time <= m->print_time + m->move_t)
            return 0.;
        m = list_next_entry(m, node);
    }
}

// Report if the given stepper is registered for the given axis
int32_t __visible
itersolve_is_active_axis(struct stepper_kinematics *sk, char axis)
{
    if (axis < 'x' || axis > 'z')
        return 0;
    return (sk->active_flags & (AF_X << (axis - 'x'))) != 0;
}

void __visible
itersolve_set_trapq(struct stepper_kinematics *sk, struct trapq *tq)
{
    sk->tq = tq;
}

void __visible
itersolve_set_stepcompress(struct stepper_kinematics *sk
                           , struct stepcompress *sc, double step_dist)
{
    sk->sc = sc;
    sk->step_dist = step_dist;
}

double __visible
itersolve_calc_position_from_coord(struct stepper_kinematics *sk
                                   , double x, double y, double z)
{
    struct move m;
    memset(&m, 0, sizeof(m));
    m.start_pos.x = x;
    m.start_pos.y = y;
    m.start_pos.z = z;
    m.move_t = 1000.;
    return sk->calc_position_cb(sk, &m, 500.);
}

void __visible
itersolve_set_position(struct stepper_kinematics *sk
                       , double x, double y, double z)
{
  double rez = itersolve_calc_position_from_coord(sk, x, y, z);
#ifdef DEBUG
  printf("Setting commanded_pos to %g\n", rez);
#endif
  sk->commanded_pos = rez;
}

double __visible
itersolve_get_commanded_pos(struct stepper_kinematics *sk)
{
    return sk->commanded_pos;
}
