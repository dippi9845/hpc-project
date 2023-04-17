#include "quad-three.h"
#include "../particle.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>

#define SQRT2 sqrtf(2.f)/2.f;

Point *newPoint(float x, float y) {
    Point * rtr = (Point *)malloc(sizeof(Point));
    rtr->x = x; rtr->y = y;
    return rtr;
}

void updateContainerForParticle(unsigned int indexParticle) {
    Point * parPos = pointFromParticle(particles + indexParticle);
    Container * myContainer = link[indexParticle];
    unsigned int containerIndex;

    for (containerIndex = 0; containerIndex < CONTAINER_CAPACITY; containerIndex++) {
        if (myContainer->particles[containerIndex] == indexParticle)
            break;
    }

    if (!isPointContained(parPos, myContainer)) {
        QuadThreeNode * validFather = myContainer->owner_node->father;

        while (validFather != NULL && !isPointContained(parPos, &(validFather->square_container))) {
            validFather = validFather->father;
        }
        
        assert(validFather != NULL); // edge case the head doesn't contains the particles ??

        while (validFather->childrens[0] != NULL)
        {
            for (int i = 0; i < CHILDREN_NUM; i++) {
                QuadThreeNode * current = validFather->childrens[i];

                if (isPointContained(parPos, &(current->square_container))) {
                    validFather = current;
                }
            }
        }
        
        // una volta trovato il nodo migliore sposta la particella li
        moveParticle(containerIndex, myContainer, &(validFather->square_container));
        
    }
    free(parPos);
}

void insertIntoContainer(Container *square, unsigned int particleIndex) {
    short int insered = 0; 

    for (int i = 0; i < CONTAINER_CAPACITY; i++) {
        if (square->particles[i] == EMPTY_INDEX) {
            insered = 1;
            square->particles[i] = particleIndex;
            link[particleIndex] = square;
            break;
        }
    }

    if (!insered) {
        // split
        Point * parPos = pointFromParticle(particles + particleIndex);
        QuadThreeNode * owner = square->owner_node;
        splitQuadNode(owner);

        for (int j = 0; j < CHILDREN_NUM; j++) {
            QuadThreeNode * current = owner->childrens[j];

            if (isPointContained(parPos, &(current->square_container))) {
                insertIntoContainer(&(current->square_container), particleIndex);
            }
        }
        free(parPos);
    }
}

void moveParticle(unsigned int particleIndex, Container *from, Container *to) {
    unsigned int particle = from->particles[particleIndex];
    from->particles[particleIndex] = EMPTY_INDEX;
    link[particleIndex] = to;
    insertIntoContainer(to, particle);
}

int isPointContained(const Point * point, const Container * square) {
    return fabsf(point->x - square->center.x) <= square->l / 2 && fabsf(point->y - square->center.y) <= square->l / 2;
}

Container *newEmptyContainerByPoints(const Point points[], QuadThreeNode *owner) {
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

Container *newEmptyContainerBySide(float side, const Point *center, QuadThreeNode *owner) {
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

QuadThreeNode * newQuadNodeByPoints(QuadThreeNode * father, const Point points[]) {
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

QuadThreeNode * newQuadNodeBySide(QuadThreeNode * father, float side, const Point * center) {
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
            Point * parPos = pointFromParticle(particles + particleIndex);

            for (int j = 0; j < CHILDREN_NUM; j++) {
                QuadThreeNode * current = node->childrens[j];

                if (isPointContained(parPos, &(current->square_container))) {
                    moveParticle(i, from, &(current->square_container));
                }
            }

            free(parPos);
        }

    }
}

Container * findContainerByPoint(QuadThreeNode * head, const Point * point) {
    
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

void insertParticle(QuadThreeNode * head, const Point * position, unsigned int indexParticle) {
    Container * to_add = findContainerByPoint(head, position);
    if (to_add != NULL) {
        insertIntoContainer(to_add, indexParticle);
    }
}

Point * pointFromParticle(const particle_t * particle) {
    return newPoint(particle->x, particle->y);
}

void applyToLeafInRange(const QuadThreeNode * head ,float radius, particle_t * pivot, void (* toApply)(particle_t *, particle_t *)) {
    const float side = head->square_container.l;
    const float halfDiagonal = side * SQRT2;
    const float maxDistance = halfDiagonal + radius;
    
    for (int i = 0; i < CHILDREN_NUM; i++) {
        const QuadThreeNode * currentNode = head->childrens[i];
        const Container * currentContainer = &(currentNode->square_container);

        const float deltaX = pivot->x - currentContainer->center.x;
        const float deltaY = pivot->y - currentContainer->center.y;
        
        if (hypotf(deltaX, deltaY) < maxDistance) {
            if (currentNode->childrens[0] == NULL) {

                // apply function to all childrens
                for (int j = 0; j < CHILDREN_NUM; j++) {
                    const unsigned int index = currentContainer->particles[j]; // index of the particle
                    toApply(pivot, particles + index);
                }
            }
            else {
                // search for leaf
                applyToLeafInRange(currentNode, radius, pivot, toApply);
            }
        }
    }

}