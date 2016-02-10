#include "common.h"
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>

#define DSP_W 1000
#define DSP_H ((int)(DSP_W*(box_h/box_w)))
#define SHRINK (0.70)

void color_background(cairo_t*, int, double, double);
void draw_particles(cairo_t*, long, double, double);
void draw_links(cairo_t*, long);
void draw_box(cairo_t*, double, double);
void draw_data(cairo_t*, double, double, double);
void draw_disk(cairo_t*, long, double, double, double, double, int);
void draw_line(cairo_t*, double, double, double, double);

void graphics(int relaxing) {
    double box_w = global.box_w;
    double box_h = global.box_h;
    long nParticles = global.nParticles;
    double time = global.time;

    cairo_surface_t *surface;
    cairo_t *cr;
    double scale = SHRINK*fmin(DSP_W/box_w, DSP_H/box_h);
    char out_file[50];
    static long framenumber = 0;
    cairo_matrix_t font_matrix;

    framenumber++;
    sprintf(out_file, "Frame%06ld.svg", framenumber);
    surface = cairo_svg_surface_create(out_file, DSP_W, DSP_H);
    cr = cairo_create(surface);
    cairo_translate(cr, 0, DSP_H);
    cairo_scale( cr, 1.0, -1.0 ); //flip y axis
    cairo_translate(cr, (1-SHRINK)/2*DSP_W, (1-SHRINK)/2*DSP_H);
    cairo_scale(cr, scale, scale);
    cairo_set_font_size(cr, fmin(box_h, box_w)*0.02);
    cairo_get_font_matrix(cr, &font_matrix);
    font_matrix.yy = -font_matrix.yy; // negative size to vertically flip text
    cairo_set_font_matrix(cr, &font_matrix);

    draw_data(cr, box_w, box_h, time);
    color_background(cr, relaxing, box_w, box_h);
    draw_particles(cr, nParticles, box_w, box_h);
    //draw_box(cr, box_w, box_h);
    draw_links(cr, nParticles);

    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    return;
}

void draw_particles(cairo_t *cr, long nParticles, double box_w,
                    double box_h) {
    cairo_set_line_width(cr,fmin(box_w, box_h)/800);
    for (long i = 0; i < nParticles; i++) {
        if (particle[i].type == 0 || particle[i].type == 1) {
            draw_disk(cr, i, particle[i].x0, particle[i].y0,
                      particle[i].radius, particle[i].w0, particle[i].type);
            draw_disk(cr, i, particle[i].x0 + box_w, particle[i].y0,
                      particle[i].radius, particle[i].w0, particle[i].type);
            draw_disk(cr, i, particle[i].x0 - box_w, particle[i].y0,
                      particle[i].radius, particle[i].w0, particle[i].type);
        } else if (particle[i].type == 2) {
            draw_line(cr, particle[i].p1x, particle[i].p1y,
                      particle[i].p2x, particle[i].p2y);
        }
    }
    return;
}

void draw_disk(cairo_t *cr, long id, double x, double y, double r,
               double angle, int type) {

    double maxrot = angle;
    double intensity = angle/maxrot;
    char label[50];
    intensity = 0.5;
    /*draw and color the disk*/
    cairo_arc(cr, x, y, r, 0, 2*M_PI);

    if (type != 0) {
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_fill(cr);
        return;
    }

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_fill_preserve(cr);

    if (angle < 0) {
        cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, intensity);
        cairo_fill_preserve(cr);
    } else if (angle > 0) {
        cairo_set_source_rgba(cr, 0.0, 0.0, 1.0, intensity);
        cairo_fill_preserve(cr);
    }

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_stroke(cr);

    /*draw a line to indicate rotation*/
    cairo_arc(cr, x, y, r, angle, angle);
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_line_to(cr, x, y);
    cairo_stroke(cr);

    //draw label
    sprintf(label,"%-3ld",id+1);
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, label);
    cairo_new_path(cr);

    return;
}

void draw_line(cairo_t *cr, double p1x, double p1y, double p2x, double p2y) {
    double box_h = global.box_h;
    double box_w = global.box_w;
    double m = (p2y - p1y)/(p2x - p1x);
    double b = p1y - m*p1x;

    if (m == 0) {
       cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
       cairo_move_to(cr, -2*box_w, b);
       cairo_line_to(cr, 2*box_w, b);
       cairo_stroke(cr);
    } else if (isinf(m)) {
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_move_to(cr, p1x, -2*box_h);
        cairo_line_to(cr, p1x, 2*box_h);
        cairo_stroke(cr);
    } else {
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_move_to(cr, -b/m, 0);
        cairo_line_to(cr, (box_h - b)/m, box_h);
        cairo_stroke(cr);
    }
    return;
}

void color_background(cairo_t *cr, int relaxing, double box_w, double box_h) {

    if (relaxing) {
        cairo_rectangle(cr, -box_w, 0, 3*box_w, box_h);
        cairo_set_source_rgb(cr, 100.0/255, 200.0/255, 100.0/255);
        cairo_fill(cr);
    }

    return;
}

void draw_box(cairo_t *cr, double box_w, double box_h) {
    cairo_set_line_width(cr,fmin(box_w, box_h)/400/2);
    cairo_rectangle(cr, 0, 0, box_w, box_h);
    cairo_set_source_rgb(cr, 1, 0, 0);
    cairo_stroke(cr);

    return;
}

void draw_data(cairo_t *cr, double box_w, double box_h, double time) {

    //draw gravity arrow
    double arrow_x = box_w + (1-SHRINK)/4*box_w/SHRINK;
    double arrow_y =  -(1-SHRINK)/4*box_h/SHRINK;
    double arrow_r =  fmin((1-SHRINK)/4*box_w/SHRINK, -arrow_y)*0.95;
    double arrow_angle = -M_PI/2-global.gravityAngle;
    double head_angle = M_PI/3.5 + M_PI - arrow_angle;
    double head_angle2 = -M_PI/3.5 + 2*M_PI - arrow_angle;
    double head_r = arrow_r*0.25;
    char line[1000];

    cairo_set_source_rgb(cr, 128/255.0, 128/255.0, 128/255.0);
    cairo_set_line_join (cr, CAIRO_LINE_JOIN_ROUND);
    cairo_set_line_width(cr,fmin(box_w, box_h)/100);
    cairo_arc(cr, arrow_x, arrow_y, arrow_r, 0, 2*M_PI);
    cairo_fill(cr);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_move_to(cr, arrow_x, arrow_y);
    cairo_arc(cr, arrow_x, arrow_y, arrow_r ,arrow_angle, arrow_angle);
    cairo_rel_line_to(cr, head_r*sin(head_angle2), head_r*cos(head_angle2));
    cairo_new_sub_path(cr);
    cairo_arc(cr, arrow_x, arrow_y, arrow_r ,arrow_angle, arrow_angle);
    cairo_rel_line_to(cr, head_r*sin(head_angle),head_r*cos(head_angle));
    cairo_stroke(cr);
    cairo_set_line_join (cr, CAIRO_LINE_JOIN_MITER);

    //Write parameters
    cairo_move_to(cr, 0, -(1-SHRINK)/8*box_h/SHRINK);
    sprintf(line, "Parámetros: Tiempo = %lf", time);
    //cairo_show_text(cr, line);
    cairo_new_path(cr);

    return;
}

void draw_links(cairo_t *cr, long nParticles) {
    long i, j;

    for (i = 0; i < nParticles; i++) {
        for (j = i + 1; j < nParticles; j++) {
            long l = i*global.nParticles + (j-1) - i*(i+3)/2;
            if (nStats[l].touching == 0) continue;
            if (particle[i].type + particle[j].type <= 1) {
                cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
                cairo_move_to(cr, particle[i].x0, particle[i].y0);
                cairo_line_to(cr, particle[j].x0, particle[j].y0);
                cairo_stroke(cr);
            } else if (particle[i].type + particle[j].type == 2
                       && particle[i].type != 1) {
                /* Wall index > disk index. Find the closest point on
                   the line (p1x,p1y)-(p2x,p2y) using the projected
                   length and subtract.*/
                double lDx = particle[j].p2x - particle[j].p1x;
                double lDy = particle[j].p2y - particle[j].p1y;
                double lenLineSqrd = lDx*lDx + lDy*lDy;
                double proy = ((particle[i].x0 - particle[j].p1x)*lDx +
                               (particle[i].y0 - particle[j].p1y)*lDy)
                    / (lenLineSqrd);
                double p2x = particle[j].p1x + proy*(lDx);
                double p2y = particle[j].p1y + proy*(lDy);
                cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
                cairo_move_to(cr, particle[i].x0, particle[i].y0);
                cairo_line_to(cr, p2x, p2y);
                cairo_stroke(cr);
            }
        }
    }
}
