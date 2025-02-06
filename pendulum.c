/*
 * Single/double pendulum simulation using X11
 *
 * Copyright 1991-2025 Chris Johns
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>
#include <X11/Xaw/Cardinals.h>
#include <X11/Xaw/XawInit.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Paned.h>
#include <X11/Xaw/Box.h>
#include <X11/Xaw/Paned.h>
#include <X11/cursorfont.h>

#define PEND1_DATA "./pend1.data"
#define PEND2_DATA "./pend2.data"

#define MAXV 50.0    /* max allowed angular velocity (radians per second) */
#define BOXLEN 300   /* length of half a box edge = radius of biggest circle */
#define PHASEBOXLEN 200 /* length of side of a phase space box */

int DEBUG = 0;
int EDEBUG = 0;

#ifdef DEBUG
#define DPRINT(fmt, args...)    printf(fmt, ## args)
#else
#define DPRINT(fmt, args...)
#endif

float t_old_c, p_old_c; /* theta/phi at end of last computation step */
float t_old_s, p_old_s; /* theta/phi at end of last display step */
float t_new, p_new;
float tdot_old, pdot_old;
float tdot_new, pdot_new;
float r1, r2, grav, lscale, m1, m2;
float dt, shdt;
float st, ct, sp, cp, sd, cd;
float dt_min = 0.003;
float dt_max = 0.08;
float initial_energy = 0.0;
int prev_time, current_time;

int boxlen, radius;
int x1_old, y1_old, x2_old, y2_old;
unsigned long foreground_pixel, background_pixel, border_pixel;
struct tms dummy_tm;
int phaseboxlen;

char *progname;
Widget box, drawform;
XtIntervalId ivid;
Widget toplevel, pane, ibox, quitbutton, resetbutton;
Widget bbox, drawbutton, stopbutton, clearbutton;
Widget tracebutton, phasebutton, stepbutton, typebutton;
Widget simbutton;
Widget phasepane1, phasepane2, phasebox1, phasebox2, phasedrawform1, phasedrawform2;
XtAppContext context;
GC draw_gc, thick_gc, xor_gc;
Window draw_win;
Display *draw_d;
Screen *draw_screen;
XtWorkProcId WorkProcId;
int screen_num;
int trace_on, phase_on, first_plot, double_flag, real_flag;
int draw_active;

extern int get_colors(void);

void
quit_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void stop_callback();

    DPRINT("quit_callback called\n");
    stop_callback(stopbutton, 0, 0);
    XtDestroyWidget(toplevel);
    DPRINT("quit_callback ended\n");
    exit(0);
}

void
stop_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    DPRINT("stop_callback called\n");
    DPRINT("ivid %d\n", ivid);
    DPRINT("draw_active = %d, WorkProcId 0x%x\n", draw_active, WorkProcId);
    DPRINT("&ivid = 0x%x, &draw_active = 0x%x\n", &ivid, &draw_active);

    if(ivid) {
        XtRemoveTimeOut(ivid);
    }
    if (draw_active) {
        XtRemoveWorkProc(WorkProcId);
    }
    XtSetSensitive(stopbutton, FALSE);
    XtSetSensitive(drawbutton, TRUE);
    XtSetSensitive(stepbutton, TRUE);
    draw_active = 0;
    DPRINT("stop_callback ended\n");
}

void
clear_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void show_callback();
    int tmp;

    DPRINT("clear_callback called\n");
    DPRINT("r1=%f r2=%f lscale=%f radius=%f\n", r1, r2, lscale, radius);

    XClearWindow(draw_d, draw_win);
    XDrawArc(draw_d, draw_win, thick_gc, 0, 0, 
            2 * boxlen, 2 * boxlen, 0, 360*64);
    first_plot = 1;
    show_callback(stepbutton, 0, 0);

    DPRINT("clear_callback ended\n");
}

void
reset_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void stop_callback();
    void clear_callback();
    void read_single_data();
    void read_double_data();
    int compute_results();

    DPRINT("reset_callback called\n");

    stop_callback(w, closure, call_data);
    if (double_flag)
        read_double_data();
    else
        read_single_data();
    (void) compute_results();
    clear_callback(w, closure, call_data);

    DPRINT("reset_callback ended\n");
}

void
trace_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    Arg arg[1];

    trace_on = 1 -  trace_on;
    if(trace_on) {
        XtSetArg(arg[0], XtNlabel, "on ");
        XtSetValues(w, arg, 1);
    } else {
        XtSetArg(arg[0], XtNlabel, "off");
        XtSetValues(w, arg, 1);
    }
}

void
phase_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    Arg arg[1];

    phase_on = 1 -  phase_on;
    if(phase_on) {
        XtSetArg(arg[0], XtNlabel, "on ");
        XtSetValues(w, arg, 1);
    } else {
        XtSetArg(arg[0], XtNlabel, "off");
        XtSetValues(w, arg, 1);
    }
}


void
real_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    Arg arg[1];
    void double_callback();

    DPRINT("real_callback called\n");

    real_flag = 1 - real_flag;

    if(real_flag) {
        XtSetArg(arg[0], XtNlabel, "real");
        XtSetValues(w, arg, 1);
    } else {
        XtSetArg(arg[0], XtNlabel, "fake");
        XtSetValues(w, arg, 1);
    }
    double_callback(typebutton, 0, 0);

    DPRINT("real_callback ended\n");
}


void
single_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void reset_callback();

    DPRINT("single_callback called\n");

    double_flag = 0;
    first_plot = 1;

    XtSetSensitive(drawbutton, TRUE);
    XtSetSensitive(stepbutton, TRUE);
    XtSetSensitive(simbutton, FALSE);

    reset_callback(resetbutton, 0, 0);

    DPRINT("single_callback ended\n");
}
    
void
type_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    Arg arg[1];
    void single_callback();
    void double_callback();

    DPRINT("type_callback called\n");

    double_flag = 1 - double_flag;
    if(double_flag) {
        XtSetArg(arg[0], XtNlabel, "double");
        XtSetValues(w, arg, 1);
        double_callback(w, closure, call_data);
    } else {
        XtSetArg(arg[0], XtNlabel, "single");
        XtSetValues(w, arg, 1);
        single_callback(w, closure, call_data);
    }

    DPRINT("type_callback ended\n");
}

void
double_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void reset_callback();

    DPRINT("double_callback called\n");

    double_flag = 1;
    first_plot = 1;

    XtSetSensitive(drawbutton, TRUE);
    XtSetSensitive(stepbutton, TRUE);
    XtSetSensitive(typebutton, TRUE);
    XtSetSensitive(simbutton, TRUE);

    reset_callback(resetbutton, 0, 0);

    DPRINT("double_callback ended\n");
}


void
draw_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void show_callback();
    int compute_results();

    DPRINT("draw_callback called\n");

    if (!draw_active) {
        XtSetSensitive(stopbutton, TRUE);
        XtSetSensitive(drawbutton, FALSE);
        XtSetSensitive(stepbutton, FALSE);
        draw_active = 1;
        prev_time = times(&dummy_tm);
    }
    WorkProcId = XtAppAddWorkProc(context, (XtWorkProc) compute_results, 0);
    show_callback(stepbutton, 0, 0);
    ivid = XtAppAddTimeOut(context, (int)(shdt*1000), show_callback, 0);

    DPRINT("draw_callback ended\n");
}

void
step_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    void show_callback();
    int compute_results();

    DPRINT("step_callback called\n");

    XtSetSensitive(stopbutton, FALSE);
    XtSetSensitive(drawbutton, TRUE);
    XtSetSensitive(stepbutton, TRUE);
    dt = shdt;  /* pretend the timestep is the time between shows */
    (void) compute_results();
    show_callback(stepbutton, 0, 0);

    DPRINT("step_callback ended\n");
}

int
compute_results(void)
{
    float tdotdot, pdotdot, A, B, C, D, E, F;
    float mu;

    if (draw_active) {
        current_time = times(&dummy_tm);
        DPRINT("compute: draw_active: current-prev = %d\n",
                current_time - prev_time);
        if ((current_time - prev_time) == 0) {
            return 0;
        }
        dt = (float) (current_time - prev_time) / 1000.0;   /* seconds */
    }
    st = sin(t_old_c);
    ct = cos(t_old_c);
    sp = sin(p_old_c);
    cp = cos(p_old_c);
    sd = sin(t_old_c - p_old_c);
    cd = cos(t_old_c - p_old_c);

    if (!double_flag) {
        tdot_new = tdot_old - dt * grav / r1 * st;
        t_new = t_old_c + tdot_new * dt;
        t_old_c = t_new;
        tdot_old = tdot_new;
    } else {
        if(!real_flag) {
            tdot_new = tdot_old - dt * grav / r1 * st;
            t_new = t_old_c + tdot_new * dt;
            pdot_new = pdot_old - dt * grav / r2 * sp;
            p_new = p_old_c + pdot_new * dt;
        } else {
            sd = sin(p_old_c - t_old_c);
            cd = cos(p_old_c - t_old_c);
            mu = m2 / m1;
            A = grav * st - mu * grav * cp * sd;
            B = mu * sd * (r2 * pdot_old * pdot_old +
                r1 * tdot_old * tdot_old * cd);
            C = (1 + mu * sd * sd) * r2;
            D = (1 + mu) * (grav * ct * sd + 
                r1 * tdot_old * tdot_old * sd);
            E = mu * r2 * pdot_old * pdot_old * sd * cd;
            F = (1 + mu * sd * sd) * r2;

            tdotdot = (B - A) / C;
            pdotdot = -(D + E) / F;

            tdot_new = tdot_old + dt * tdotdot;
            pdot_new = pdot_old + dt * pdotdot;

            t_new = t_old_c + dt * tdot_new;
            p_new = p_old_c + dt * pdot_new;
        }
        t_old_c = t_new;
        tdot_old = tdot_new;
        p_old_c = p_new;
        pdot_old = pdot_new;
    }
    prev_time = current_time;
    return 0;
}

void
show_callback(Widget w, caddr_t closure, caddr_t call_data)
{
    int x1_new, y1_new, x2_new, y2_new;
    float kinetic, potential, energy, delta_e;

    st = sin(t_old_s);
    ct = cos(t_old_s);
    sp = sin(p_old_s);
    cp = cos(p_old_s);
    sd = sin(t_old_s - p_old_s);
    cd = cos(t_old_s - p_old_s);

    DPRINT("show_callback called\n");

    if (first_plot) {
        x1_old = (int) (boxlen + lscale * r1 * st);
        y1_old = (int) (boxlen + lscale * r1 * ct);
        x2_old = (int) (x1_old + lscale * r2 * sp);
        y2_old = (int) (y1_old + lscale * r2 * cp);
    }

    if (!double_flag) {

        kinetic = 0.5 * r1 * r1 * tdot_old * tdot_old;
        potential = grav * r1 * (1.0 - ct);

        DPRINT("kinetic = %15.5f, potential = %15.5f, total = %15.5f\n",
                kinetic, potential, kinetic + potential);

        if(first_plot) {
            XDrawLine(draw_d, draw_win, xor_gc, boxlen, boxlen, x1_old, y1_old);
            XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x1_old - radius),
                (int) (y1_old - radius),
                2 * radius, 2 * radius, 0, 360*64);
            if(trace_on) {
                XDrawPoint(draw_d, draw_win, draw_gc, 
                (int) x1_old, (int) y1_old);
            }
            first_plot = 0;
        }

        XDrawLine(draw_d, draw_win, xor_gc, boxlen, boxlen, x1_old, y1_old);
        XFillArc(draw_d, draw_win, xor_gc, 
            (int) (x1_old - radius),
            (int) (y1_old - radius),
            2 * radius, 2 * radius, 0, 360*64);

        x1_new = (int) (boxlen + lscale * r1 * sin(t_new));
        y1_new = (int) (boxlen + lscale * r1 * cos(t_new));

        if (trace_on && !first_plot) {
            XDrawLine(draw_d, draw_win, draw_gc, 
                (int) x1_old, (int) y1_old,
                (int) x1_new, (int) y1_new);
        }
        x1_old = x1_new;
        y1_old = y1_new;
        XDrawLine(draw_d, draw_win, xor_gc, boxlen, boxlen, x1_old, y1_old);
        XFillArc(draw_d, draw_win, xor_gc, 
            (int) (x1_old - radius),
            (int) (y1_old - radius),
            2 * radius, 2 * radius, 0, 360*64);
    } else {

        if (!real_flag) {
            kinetic = 0.5 * r1 * r1 * tdot_old * tdot_old
                    + 0.5 * r2 * r2 * pdot_old * pdot_old;
            potential = grav * r1 * (1.0 - ct)
                    + grav * r2 * (1.0 - cp);
            DPRINT("kinetic = %15.5f, potential = %15.5f, total = %15.5f\n",
                    kinetic, potential, kinetic + potential);
        } else {
            kinetic = 0.5 * m1 * r1 * r1 * tdot_old * tdot_old +
                0.5 * m2 * (r1 * r1 * tdot_old * tdot_old +
                r2 * r2 * pdot_old * pdot_old + 2 * r1 * r2 * tdot_old
                * pdot_old * cd),
            potential = m1 * grav * (r1 + r2 - r1 * ct) +
                m2 * grav * (r1 + r2 - (r1 * ct + r2 * cp));

            DPRINT("kinetic = %15.5f, potential = %15.5f, total = %15.5f\n",
                    kinetic, potential, kinetic + potential);
        }

        if(first_plot) {
            XDrawLine(draw_d, draw_win, xor_gc, boxlen, boxlen, 
                    x1_old, y1_old);
            XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x1_old - radius),
                (int) (y1_old - radius),
                2 * radius, 2 * radius, 0, 360*64);
            XDrawLine(draw_d, draw_win, xor_gc, x1_old, y1_old, 
                    x2_old, y2_old);
            XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x2_old - radius),
                (int) (y2_old - radius),
                2 * radius, 2 * radius, 0, 360*64);
            if(trace_on) {
                XDrawPoint(draw_d, draw_win, draw_gc, 
                    (int) x2_old, (int) y2_old);
            }
            first_plot = 0;
        }
        XDrawLine(draw_d, draw_win, xor_gc, boxlen, boxlen, x1_old, y1_old);
        XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x1_old - radius),
                (int) (y1_old - radius),
                2 * radius, 2 * radius, 0, 360*64);
        XDrawLine(draw_d, draw_win, xor_gc, x1_old, y1_old, x2_old, y2_old);
        XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x2_old - radius),
                (int) (y2_old - radius),
                2 * radius, 2 * radius, 0, 360*64);


        x1_new = (int) (boxlen + lscale * r1 * st);
        y1_new = (int) (boxlen + lscale * r1 * ct);
        x2_new = (int) (x1_new + lscale * r2 * sp);
        y2_new = (int) (y1_new + lscale * r2 * cp);
        if (trace_on && !first_plot) {
            XDrawLine(draw_d, draw_win, draw_gc, 
                    (int) x2_old, (int) y2_old,
                    (int) x2_new, (int) y2_new);
        }
        x1_old = x1_new;
        y1_old = y1_new;
        x2_old = x2_new;
        y2_old = y2_new;
        XDrawLine(draw_d, draw_win, xor_gc, boxlen, boxlen, x1_old, y1_old);
        XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x1_old - radius),
                (int) (y1_old - radius),
                2 * radius, 2 * radius, 0, 360*64);
        XDrawLine(draw_d, draw_win, xor_gc, x1_old, y1_old, x2_old, y2_old);
        XFillArc(draw_d, draw_win, xor_gc, 
                (int) (x2_old - radius),
                (int) (y2_old - radius),
                2 * radius, 2 * radius, 0, 360*64);
    }
    t_old_s = t_new;
    p_old_s = p_new;

    DPRINT("show_callback ended\n");

    if (draw_active)
        ivid = XtAppAddTimeOut(context, (int)(shdt*1000),
                    (XtTimerCallbackProc)show_callback, 0);
}


void
read_single_data()
{
    FILE *datafile;
    char buffer[80];

    if ((datafile = fopen(PEND1_DATA, "r")) == NULL) {
        fprintf(stderr, "Unable to open data file %s\n", PEND1_DATA);
        exit(1);
    }
    fscanf(datafile, "%s\n", buffer);
    t_old_c = t_old_s = atof(buffer);

    fscanf(datafile, "%s\n", buffer);
    tdot_old = atof(buffer);

    fscanf(datafile, "%s\n", buffer);
    grav = atof(buffer);

    fscanf(datafile, "%s\n", buffer);
    r1 = atof(buffer);

    fscanf(datafile, "%s\n", buffer);
    shdt = atof(buffer);

    fscanf(datafile, "%s\n", buffer);
    radius = atoi(buffer);

    fclose(datafile);

    r2 = m1 = m2 = 0.0;
    boxlen = BOXLEN;
    lscale = (boxlen - radius) / (r1 + r2);
    phaseboxlen = PHASEBOXLEN;
}


void
read_double_data()
{
    FILE *datafile;
    char buffer[80];

    if ((datafile = fopen(PEND2_DATA, "r")) == NULL) {
        fprintf(stderr, "Unable to open data file %s\n", PEND2_DATA);
        exit(1);
    }
    fscanf(datafile, "%s\n", buffer);
    t_old_c = t_old_s = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    tdot_old = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    p_old_c = p_old_s = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    pdot_old = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    grav = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    r1 = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    r2 = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    shdt = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    m1 = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    m2 = (float) strtod(buffer, (char **)0);

    fscanf(datafile, "%s\n", buffer);
    radius = (float) strtod(buffer, (char **)0);

    fclose(datafile);

    boxlen = BOXLEN;
    lscale = (boxlen - radius) / (r1 + r2);
    phaseboxlen = PHASEBOXLEN;
}

void
app_warning(void)
{
    fprintf(stderr, "Warning from application!\n");
}

void
app_error(void)
{
    fprintf(stderr, "ERROR from application!\n");
}

int
main(int argc, char **argv)
{
    int i;
    Arg argies[10];
    Widget gaplabel, tracelabel, phaselabel, singlelabel, simlabel;
    int ret;

    trace_on = 0;
    phase_on = 0;
    first_plot = 1;
    double_flag = 1;
    real_flag = 1;
    draw_active = 0;

    if (double_flag)
        read_double_data();
    else
        read_single_data();

    progname = argv[0];
    toplevel = XtInitialize(argv[0], "Pendulum", NULL, 0, &argc, argv);
    draw_d = XtDisplay(toplevel);

    /* Main window */
    pane = XtCreateManagedWidget("pane", panedWidgetClass, toplevel,
        NULL, 0);

    box = XtCreateManagedWidget("box", boxWidgetClass, pane,
        NULL, 0);

    i = 0;
    XtSetArg(argies[i], XtNheight, (XtArgVal) (2 * boxlen)); i++;
    XtSetArg(argies[i], XtNwidth, (XtArgVal) (2 * boxlen)); i++;
    XtSetArg(argies[i], XtNtop, XtChainTop); i++;
    XtSetArg(argies[i], XtNbottom, XtChainBottom); i++;
    XtSetArg(argies[i], XtNleft, XtChainLeft); i++;
    XtSetArg(argies[i], XtNright, XtChainRight); i++;
    XtSetArg(argies[i], XtNforeground, (XtArgVal)
        WhitePixel(draw_d, DefaultScreen(draw_d))); i++;
    XtSetArg(argies[i], XtNbackground, (XtArgVal)
        BlackPixel(draw_d, DefaultScreen(draw_d))); i++;
    drawform = XtCreateManagedWidget("drawform", coreWidgetClass,
        pane, argies, i);

    /* Phase space window #1 */
    phasepane1 = XtCreateManagedWidget("phasepane1", panedWidgetClass, toplevel,
        NULL, 0);

    phasebox1 = XtCreateManagedWidget("phasebox1", boxWidgetClass, phasepane1,
        NULL, 0);

    i = 0;
    XtSetArg(argies[i], XtNheight, (XtArgVal) phaseboxlen); i++;
    XtSetArg(argies[i], XtNwidth, (XtArgVal) phaseboxlen); i++;
    XtSetArg(argies[i], XtNtop, XtChainTop); i++;
    XtSetArg(argies[i], XtNbottom, XtChainBottom); i++;
    XtSetArg(argies[i], XtNleft, XtChainLeft); i++;
    XtSetArg(argies[i], XtNright, XtChainRight); i++;
    XtSetArg(argies[i], XtNforeground, (XtArgVal)
        WhitePixel(draw_d, DefaultScreen(draw_d))); i++;
    XtSetArg(argies[i], XtNbackground, (XtArgVal)
        BlackPixel(draw_d, DefaultScreen(draw_d))); i++;
    phasedrawform1 = XtCreateManagedWidget("phasedrawform1", coreWidgetClass,
        phasepane1, argies, i);

    /* Phase space window #2 */
    phasepane2 = XtCreateManagedWidget("phasepane2", panedWidgetClass, toplevel,
        NULL, 0);

    phasebox2 = XtCreateManagedWidget("phasebox2", boxWidgetClass, phasepane2,
        NULL, 0);

    i = 0;
    XtSetArg(argies[i], XtNheight, (XtArgVal) phaseboxlen); i++;
    XtSetArg(argies[i], XtNwidth, (XtArgVal) phaseboxlen); i++;
    XtSetArg(argies[i], XtNtop, XtChainTop); i++;
    XtSetArg(argies[i], XtNbottom, XtChainBottom); i++;
    XtSetArg(argies[i], XtNleft, XtChainLeft); i++;
    XtSetArg(argies[i], XtNright, XtChainRight); i++;
    XtSetArg(argies[i], XtNforeground, (XtArgVal)
        WhitePixel(draw_d, DefaultScreen(draw_d))); i++;
    XtSetArg(argies[i], XtNbackground, (XtArgVal)
        BlackPixel(draw_d, DefaultScreen(draw_d))); i++;
    phasedrawform2 = XtCreateManagedWidget("phasedrawform2", coreWidgetClass,
        phasepane2, argies, i);


    i = 0;
    quitbutton = XtCreateManagedWidget("quit", commandWidgetClass,
        box, argies, i);
    XtAddCallback(quitbutton, XtNcallback,
            (XtCallbackProc)quit_callback, NULL);

    i = 0;
    XtSetArg(argies[i], XtNborderWidth, (XtArgVal) 0); i++;
    gaplabel = XtCreateManagedWidget("", labelWidgetClass,
        box, argies, i);

    i = 0;
    clearbutton = XtCreateManagedWidget("clear", commandWidgetClass,
        box, argies, i);
    XtAddCallback(clearbutton, XtNcallback,
            (XtCallbackProc)clear_callback, NULL);
    i = 0;
    stopbutton = XtCreateManagedWidget("stop", commandWidgetClass,
        box, argies, i);
    XtAddCallback(stopbutton, XtNcallback,
            (XtCallbackProc)stop_callback, NULL);
    i = 0;
    stepbutton = XtCreateManagedWidget("step", commandWidgetClass,
        box, argies, i);
    XtAddCallback(stepbutton, XtNcallback,
            (XtCallbackProc)step_callback, NULL);
    i = 0;
    drawbutton = XtCreateManagedWidget("draw", commandWidgetClass,
        box, argies, i);
    XtAddCallback(drawbutton, XtNcallback,
            (XtCallbackProc)draw_callback, NULL);
    i = 0;
    resetbutton = XtCreateManagedWidget("reset", commandWidgetClass,
        box, argies, i);
    XtAddCallback(resetbutton, XtNcallback,
            (XtCallbackProc)reset_callback, NULL);

    i = 0;
    XtSetArg(argies[i], XtNborderWidth, (XtArgVal) 0); i++;
    tracelabel = XtCreateManagedWidget("trace:", labelWidgetClass,
        box, argies, i);
    i = 0;
    tracebutton = XtCreateManagedWidget("off", commandWidgetClass,
        box, argies, i);
    XtAddCallback(tracebutton, XtNcallback,
            (XtCallbackProc)trace_callback, NULL);

    i = 0;
    XtSetArg(argies[i], XtNborderWidth, (XtArgVal) 0); i++;
    phaselabel = XtCreateManagedWidget("phase:", labelWidgetClass,
        box, argies, i);
    i = 0;
    phasebutton = XtCreateManagedWidget("off", commandWidgetClass,
        box, argies, i);
    XtAddCallback(phasebutton, XtNcallback,
            (XtCallbackProc)phase_callback, NULL);

    i = 0;
    XtSetArg(argies[i], XtNborderWidth, (XtArgVal) 0); i++;
    singlelabel = XtCreateManagedWidget("type:", labelWidgetClass,
        box, argies, i);
    i = 0;
    typebutton = XtCreateManagedWidget("double", commandWidgetClass,
        box, argies, i);
    XtAddCallback(typebutton, XtNcallback,
            (XtCallbackProc)type_callback, NULL);

    i = 0;
    XtSetArg(argies[i], XtNborderWidth, (XtArgVal) 0); i++;
    simlabel = XtCreateManagedWidget("simtype:", labelWidgetClass,
        box, argies, i);
    i = 0;
    simbutton = XtCreateManagedWidget("real", commandWidgetClass,
        box, argies, i);
    XtAddCallback(simbutton, XtNcallback,
            (XtCallbackProc)real_callback, NULL);

    XtAddEventHandler(drawform, ExposureMask, FALSE,
                (XtEventHandler)clear_callback, 0);

    XtRealizeWidget(toplevel);
    draw_d = XtDisplay(drawform);
    draw_win = XtWindow(drawform);
    draw_screen = XtScreen(drawform);
    screen_num = DefaultScreen(draw_d);
    context = XtWidgetToApplicationContext(drawform);

#if 0
    XtAppSetErrorHandler(context, (XtErrorHandler) app_error);
    XtAppSetWarningHandler(context, (XtErrorHandler) app_warning);
#endif

    ret = get_colors();
    draw_gc = draw_screen->default_gc;

    xor_gc = XCreateGC(draw_d, draw_win, 0, 0);
    XCopyGC(draw_d, draw_gc, -1, xor_gc);
    XSetFunction(draw_d, xor_gc, GXxor);

    thick_gc = XCreateGC(draw_d, draw_win, 0, 0);
    XCopyGC(draw_d, draw_gc, -1, thick_gc);
    XSetLineAttributes(draw_d, thick_gc, 3, LineSolid, CapButt, JoinMiter);

    XSetForeground(draw_d, xor_gc, foreground_pixel);
    XSetForeground(draw_d, draw_gc, foreground_pixel);
    XSetForeground(draw_d, thick_gc, foreground_pixel);
    XSetBackground(draw_d, xor_gc, background_pixel);
    XSetBackground(draw_d, draw_gc, background_pixel);
    XSetBackground(draw_d, thick_gc, background_pixel);

    if (double_flag)
        double_callback(typebutton, 0, 0);
    else
        single_callback(typebutton, 0, 0);

    XtMainLoop();
}
