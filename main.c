#include <stdio.h>
#include <stdlib.h>

typedef struct _particle {
    double x, y, z, mass;
} Particle;

typedef struct _particleList {
    int count;
    Particle *particles;
} ParticleList;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"

ParticleList parse(const char *path) {
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        printf("File %s could not be opened!\n", path);
        exit(1);
    }

    int count;
    fscanf(f, "%d\n", &count);

    Particle *particles = malloc(count * sizeof(Particle));
    for (int i = 0; i < count; i++) {
        fscanf(f, "%lf ", &particles[i].x);
        fscanf(f, "%lf ", &particles[i].y);
        fscanf(f, "%lf ", &particles[i].z);
        fscanf(f, "%lf\n", &particles[i].mass);
    }
    fclose(f);

    return (ParticleList) {
            .count = count,
            .particles = particles
    };
}

#pragma clang diagnostic pop

Particle centerOfGravity(ParticleList list) {
    double x = 0, y = 0, z = 0, totalMass = 0;
    for (int i = 0; i < list.count; i++) {
        Particle p = list.particles[i];
        x += p.mass * p.x;
        y += p.mass * p.y;
        z += p.mass * p.z;
        totalMass += p.mass;
    }
    x /= totalMass;
    y /= totalMass;
    z /= totalMass;

    return (Particle) {
            .x = x,
            .y = y,
            .z = z,
            .mass = totalMass
    };
}

int main() {
    ParticleList list = parse("lab01-input.txt");
    Particle center = centerOfGravity(list);
    printf("%e\t%e\t%e\t%e\n", center.x, center.y, center.z, center.mass);
    return 0;
}