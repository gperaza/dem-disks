#include "common.h"

int front_add();
void draw_frame(long i);
void place_disk_dd(long i);
void place_disk_dw(long i, disk*, disk*);
int check_overlap(disk*, long);
void check_overlap_all();

typedef struct Front {
    disk *first;
    disk *last;
    /* Active front, from a to b. */
    disk *a;
    disk *b;
} Front;

Front front;
int layer;

void pack_disks() {
    long nDisks = global.nParticles - 4;
    /* Builds a disk packing using an advancing front approach. */

    /* Build the initial front including walls. */
    /* Front is Wall_L(a) -> Wall_B(b) -> Wall_R. */
    front.a = front.first = &particle[nDisks+3];  // left wall
    front.a->next = front.b = &particle[nDisks];  // bottom wall
    front.a->prev = NULL;
    front.last = &particle[nDisks+1];  // right wall
    front.b->next = front.last;
    front.b->prev = front.a;
    front.last->next = NULL;
    front.last->prev = front.b;

    /* Add first disk manually and update front. */
    /* Front is Wall_L -> Disk[0](a) -> Wall_B(b) -> Wall_R. */
    particle[0].x0 = particle[0].radius;
    particle[0].y0 = particle[0].radius;
    particle[0].prev = front.a;
    particle[0].next = front.b;
    front.a->next = front.b->prev = &particle[0];
    front.a = &particle[0];

    /* Build layers. */
    int success;
    for (long i = 1; i < nDisks; i++) {
        success = 0;
        while (!success) {
            success = front_add(i);
        }
        if (front.b == front.last) {
            /* Restart front. */
            front.a = front.first;
            front.b = front.a->next;
            layer++;
        }
    }
}

int front_add(long i) {
    disk *curr;

    assert(front.a->type + front.b->type <= 2);

    /* Try to place disk. */
    if (front.a->type == 0 && front.b->type == 0) place_disk_dd(i);
    else if (front.a->type == 2) place_disk_dw(i, front.b, front.a);
    else place_disk_dw(i, front.a, front.b);


    /* Check overlap backward. */
    int overlapb = 0;
    disk *backColl = NULL;
    curr = front.a->prev;
    while (curr != NULL && !overlapb) {
        if (check_overlap(curr, i)) {
            overlapb++;
            backColl = curr;
        }
        curr = curr->prev;
    }

    /* Check overlap forward. */
    int overlapf = 0;
    disk *forwardColl = NULL;
    curr = front.b->next;
    while (curr != NULL && !overlapf) {
        if (check_overlap(curr, i)) {
            overlapf++;
            forwardColl = curr;
        }
        curr = curr->next;
    }

    /* If no overlap advance front. */
    if (overlapf + overlapb == 0) {
        front.a->next = &particle[i];
        particle[i].next = front.b;
        particle[i].prev = front.a;
        front.b->prev = &particle[i];
        if (layer == 0) front.a = &particle[i];
        else {
            front.a = front.b;
            front.b = (front.a->next != NULL) ? front.a->next : front.last;
        }
        return 1;
    }

    /* If just forward overlap adjust front. */
    if (overlapf) {
        forwardColl->prev = front.a;
        front.a->next = forwardColl;
        front.b = forwardColl;
    }

    /* If just backward overlap. */
    if (overlapb) {
        backColl->next = front.b;
        front.b->prev = backColl;
        front.a = backColl;
    }

    return 0;
}

int check_overlap(disk* curr, long i) {
    int overlap = 0;
    double dx, dy, d2, d;

    if (curr->type == 0 /*disk*/) {
        dx = particle[i].x0 - curr->x0;
        dy = particle[i].y0 - curr->y0;
        d2 = dx*dx + dy*dy;
        if (d2 < (particle[i].radius + curr->radius)*
            (particle[i].radius + curr->radius)) overlap = 1;
    } else /*wall*/ {
        dx = curr->nx*(particle[i].x0 - curr->x0);
        dy = curr->ny*(particle[i].y0 - curr->y0);
        d = dx + dy;
        if (fabs(d) < particle[i].radius) overlap = 1;
    }
    return overlap;
}

void place_disk_dd(long i) {
    /* Place new disk touching other two.*/
    double Dx = front.b->x0 - front.a->x0;
    double Dy = front.b->y0 - front.a->y0;
    double L2 = Dx*Dx + Dy*Dy;
    double Ria = particle[i].radius + front.a->radius;
    double Rib = particle[i].radius + front.b->radius;
    double A = L2 + Ria*Ria - Rib*Rib;
    double B = sqrt((L2 - (Ria - Rib)*(Ria - Rib))*
                    ((Ria + Rib)*(Ria + Rib) - L2));

    particle[i].y0 = front.a->y0 + (A*Dy + B*Dx)*0.5/L2;
    particle[i].x0 = front.a->x0 + (A*Dx - B*Dy)*0.5/L2;
}

void place_disk_dw(long i, disk *d, disk *w) {
    /* Place new disk touching a wall and a disk.*/
    /* Signs in this function take into account the fact that walls
       are orthogonal and we always want the grater solution. It is
       not a general solution.*/

    double Dx = w->x0 - d->x0;
    double Dy = w->y0 - d->y0;
    double Ria = particle[i].radius + d->radius;
    double A = w->nx*Dx + w->ny*Dy + particle[i].radius;
    double B = sqrt(Ria*Ria - A*A);

    particle[i].y0 = d->y0 + A*w->ny + B*fabs(w->nx);
    particle[i].x0 = d->x0 + A*w->nx + B*w->ny;
}
