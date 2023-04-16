#include "quad-three.h"
#include <stdlib.h>
#include <math.h>

Point *newPoint(float x, float y) {
    Point * rtr = (Point *)malloc(sizeof(Point));
    rtr->x = x; rtr->y = y;
    return rtr;
}

Container *findNearContainerForPoint(const Point *point, const Container *square) {
    return NULL;
}

void insertIntoContainer(Container *square, unsigned int particleIndex) {
    short int insered = 0; 

    for (int i = 0; i < CONTAINER_CAPACITY; i++) {
        if (square->particles[i] == EMPTY_INDEX) {
            insered = 1;
            square->particles[i] = particleIndex;
            break;
        }
    }

    if (!insered) {
        // split
        splitQuadNode(square->owner_node);
    }
}

void moveParticle(unsigned int particleIndex, Container *from, Container *to) {
    unsigned int particle = from->particles[particleIndex];
    from->particles[particleIndex] = EMPTY_INDEX;
    insertIntoContainer(to, particle);
}

int isPointContained(const Point * point, const Container * square) {
    return fabsf(point->x - square->center.x) <= square->l / 2 && fabsf(point->y - square->center.y) <= square->l / 2;
}

Container *newEmptyContainerByPoints(const Point points[], const QuadThreeNode *owner) {
    Container * rtr = (Container *) malloc(sizeof(Container));
    
    rtr->owner_node = owner;

    for (int i = 0; i < CHILDREN_NUM; i++) {
        rtr->points[i] = points[i];
    }

    rtr->l = fabsf(points[0].x - points[1].x);

    Point *tmp = newPoint(points[TOP_LEFT].x + rtr->l / 2, points[TOP_LEFT].y + rtr->l / 2);
    rtr->center = *tmp;
    free(tmp); // avoids memory leak

    for (int i = 0; i < CONTAINER_CAPACITY; i++) {
        rtr->particles[i] = EMPTY_INDEX;
    }

    return rtr;
}

Container *newEmptyContainerBySide(float side, const Point *center, const QuadThreeNode *owner) {
    Container * rtr = (Container *) malloc(sizeof(Container));
    
    rtr->owner_node = owner;
    rtr->l = side;
    rtr->center = *center;

    for (int i = 0; i < CONTAINER_CAPACITY; i++) {
        rtr->particles[i] = EMPTY_INDEX;
    }

    Point *tmp = newPoint(center->x - side/2, center->y - side/2);
    
    rtr->points[TOP_LEFT] = *tmp;
    free(tmp); // avoids memory leak

    tmp = newPoint(center->x + side/2, center->y - side/2);

    rtr->points[TOP_RIGHT] = *tmp;
    free(tmp);

    tmp = newPoint(center->x - side/2, center->y + side/2);

    rtr->points[BOTTOM_LEFT] = *tmp;
    free(tmp);

    tmp = newPoint(center->x + side/2, center->y + side/2);

    rtr->points[BOTTOM_RIGHT] = *tmp;
    free(tmp);

    return rtr;
}

QuadThreeNode * newQuadNodeByPoints(const QuadThreeNode * father, const Point points[]) {
    QuadThreeNode * rtr = (QuadThreeNode *) malloc(sizeof(QuadThreeNode));
    
    for (int i = 0; i < CHILDREN_NUM; i++) {
        rtr->childrens[i] = NULL;
    }

    rtr->father = father;
    
    Container * tmp = newEmptyContainerByPoints(points, rtr);
    
    rtr->square_container = *tmp;
    
    free(tmp); // avoids memory leak

    return rtr;
}

QuadThreeNode * newQuadNodeBySide(const QuadThreeNode * father, float side, const Point * center) {
    QuadThreeNode * rtr = (QuadThreeNode *) malloc(sizeof(QuadThreeNode));
    
    for (int i = 0; i < CHILDREN_NUM; i++) {
        rtr->childrens[i] = NULL;
    }

    rtr->father = father;
    
    Container * tmp = newEmptyContainerBySide(side, center, rtr);
    
    rtr->square_container = *tmp;
    
    free(tmp); // avoids memory leak
    
    return rtr;
}

void splitQuadNode(QuadThreeNode * node) {
    if (node->childrens[0] == NULL) {
        float side = node->square_container.l;
        Point * center = &(node->square_container.center);

        Point * tmp = newPoint(center->x - side/4, center->y - side/4);

        node->childrens[TOP_LEFT] = newQuadNodeBySide(node, side/2, tmp);
        free(tmp); // avoids memory leak

        tmp = newPoint(center->x + side/4, center->y - side/4);

        node->childrens[TOP_RIGHT] = newQuadNodeBySide(node, side/2, tmp);
        free(tmp);

        tmp = newPoint(center->x - side/4, center->y + side/4);

        node->childrens[BOTTOM_LEFT] = newQuadNodeBySide(node, side/2, tmp);
        free(tmp);

        tmp = newPoint(center->x + side/4, center->y + side/4);

        node->childrens[BOTTOM_RIGHT] = newQuadNodeBySide(node, side/2, tmp);
        free(tmp);

        // sposta le particelle
        Container * from = &(node->square_container);

        for (int i = 0; i < CONTAINER_CAPACITY; i++) {
            unsigned int particleIndex = from->particles[i];

            for (int j = 0; j < CHILDREN_NUM; j++) {
                const QuadThreeNode * current = node->childrens[j];
                Point * parPos = pointFromParticle(particles[particleIndex]);

                if (isPointContained(parPos, &(node->square_container))) {
                    moveParticle(i, from, &(current->square_container));
                }
            }
        }

    }
}

Container * findContainerByPoint(const QuadThreeNode * head, const Point * point) {
    
    if (head->childrens[0] == NULL)
        return &(head->square_container);
    
    for (int i = 0; i < CHILDREN_NUM; i++) {

        QuadThreeNode * currentNode = head->childrens[i];
        
        if (isPointContained(point, &(currentNode->square_container))) {
            return findContainerByPoint(currentNode, point);
        }
    }

    return NULL;
}

void insertParticle(const QuadThreeNode * head, const Point * position, unsigned int indexParticle) {
    Container * to_add = findContainerByPoint(head, position);
    if (to_add != NULL) {
        insertIntoContainer(to_add, indexParticle);
    }
}

Point * pointFromParticle(const particle_t * particle) {
    return newPoint(particle->x, particle->y);
}