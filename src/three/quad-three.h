#include "../particle.h"

#ifndef QUAD_THREE
#define QUAD_THREE

#define TOP_LEFT 0
#define TOP_RIGHT 1
#define BOTTOM_LEFT 2
#define BOTTOM_RIGHT 3

#define CONTAINER_CAPACITY 10
#define EMPTY_INDEX 0xFFFFFFFF

#define CHILDREN_NUM 4

typedef struct {
    float x, y;
} Point;

/* return a point from given coordinates */
Point * newPoint(float x, float y);

inline Point * newPointZero() { return newPoint( 0.f, 0.f ); }


typedef struct {
    float l; // lato
    Point center;
    Point points[4];                     // 4 spigoli del quadrato
    unsigned int particles[CONTAINER_CAPACITY]; // indici delle particelle contenute
    struct QuadThreeNode * owner_node;

} Container;

/* links the particle index to a container that contains it particle */
Container ** link;

/* checks if the current container of this particle is still valid, otherwise it will be changed form one near */
void updateContainerForParticle(unsigned int indexParticle);

/* insert the index into the container */
void insertIntoContainer(Container * square, unsigned int particleIndex);

/* particleIndex: refers to the index of the particle of the from container
    moves the particle to the container "to"
*/
void moveParticle(unsigned int particleIndex, Container * from, Container * to);

/* return True if point is in range of the container */
int isPointContained(const Point * point, const Container * square);


typedef struct QuadThreeNode {
    struct QuadThreeNode * childrens[4];
    struct QuadThreeNode * father;
    Container square_container;
} QuadThreeNode;


Container * newEmptyContainerByPoints(const Point points[],  const QuadThreeNode * owner);
Container * newEmptyContainerBySide(float side, const Point * center, const QuadThreeNode * owner);

QuadThreeNode * newQuadNodeByPoints(const QuadThreeNode * father, const Point points[]);
QuadThreeNode * newQuadNodeBySide(const QuadThreeNode * father, float side, const Point * center);

/* the childrens of the ginven node will be defined, if they aren't */
void splitQuadNode(QuadThreeNode * node);

/* find container that would */
Container * findContainerByPoint(const QuadThreeNode * head, const Point * point);

/* insert index into quadthree */
void insertParticle(const QuadThreeNode * head, const Point * position, unsigned int indexParticle);

/* returns a point thats have the positon of the particle given */
Point * pointFromParticle(const particle_t * particle);

#endif