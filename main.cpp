#include <QApplication>

#include <iostream>

#include "camera.h"
#include "gui.h"

// for test
#include "ray.h"
#include "entities.h"
#include "glm/ext.hpp"

// test functions
void entity_test();
void bbox_test();
void matrix_test();
void insert_tris(Octree& scene);

int main(int argc, char** argv) {
    QApplication app(argc, argv);
//    entity_test();

    //Camera camera({10, 0, 0});
    Camera camera({-10, 0, 0}, {1, 0, 0}, 0.1);
    glm::dvec3 light{-10, 10, 10};

    RayTracer raytracer(camera, light);

    // Set up scene
    Octree scene({-20, -20, -20}, {20, 20, 20});
//    ImpSphere *s = new ImpSphere(glm::dvec3{4,0,0}, 2, {0,1,0});
//    ImpTriangle *r = new ImpTriangle({3,-2,1},{3,2,-1}, {3,-2,-1});
    ImpSphere *s2 = new ImpSphere(glm::dvec3{3,4,4}, 2, {1,0,0});
    ImpSphere *s3 = new ImpSphere(glm::dvec3{4,-4,4}, 2, {0,0,1});
    
//    ImpTriangle *t = new ImpTriangle({0,3,2},{3,3,-4},{3,-3,-4});
//    ExpCone *c = new ExpCone({0,0,2}, {-1,1,-3}, 5, 3, {1,1,0}); //-1,1,-3
//    ExpSphere *s = new ExpSphere(glm::dvec3{-2,0,0}, 2, {0,1,0});
//    ExpCube *cube = new ExpCube(glm::dvec3{0,0,0},2,2,2,{1,0,0});
//    ExpRectangle *r = new ExpRectangle(glm::dvec3{0,0,0}, glm::dvec3{3,3,3}, glm::dvec3{3,3,0});
    
    ExpQuad *q = new ExpQuad(glm::dvec3{0,0,0},2,3, (90.0* M_PI / 180.0),{1,2,3});  // rotate around y by 80 degree


//    scene.push_back(s);
    scene.push_back(q);
    scene.push_back(s2);
    scene.push_back(s3);
    
//    scene.push_back(q);
//    scene.push_back(s);
//    insert_tris(scene);
    
    // TODO Add objects to the scene
    // scene.push_back(...);

    raytracer.setScene(&scene);

    Gui window(500, 500, raytracer);
    window.show();
    return app.exec();
}

void insert_tris(Octree& scene) {
    glm::dvec3 pos = {0,0,0};
    float alpha;
    std::vector<glm::dvec3> vertices;
    glm::dvec3 loc = glm::dvec3{0,0,0};
    vertices.push_back(pos);
    double radius = 1;
    double height = 2;
    
    double numSubdivisions = 23.0;
    for (int i = 0; i <= numSubdivisions; ++i) {
        alpha = i * 360.0/numSubdivisions;
        loc.x = pos.x + radius * cos(alpha * PI / 180.0);
        //            std::cout << "alpha: " << alpha << std::endl;
        //            std::cout << "cos(alpha): " << cos(alpha * PI / 180.0) << std::endl;
        loc.y = pos.y + radius * sin(alpha * PI / 180.0);
        loc.z = pos.z - height;
        vertices.push_back(loc);
    }
    
//    for(int i = 1; i < vertices.size()-15; ++i){
//        scene.push_back(new ImpTriangle(vertices[0],vertices[i],vertices[i+1]));
//    }
    scene.push_back(new ImpTriangle(vertices[0],vertices[1],vertices[2]));
}

void entity_test() {
    ImpSphere s = ImpSphere(glm::dvec3{2,0,0}, 10, {0,1,0});
    //std::cout << 2 / 0 << std::endl;
    std::cout << "sphere created" << std::endl;
    std::cout << "sphere radius: " << s.radius << std::endl;
    std::cout << "sphere position: " << glm::to_string(s.pos) << std::endl;
    
    Ray r = Ray(glm::dvec3{-10, 0, 0}, glm::dvec3{1,0.5,0.5});
    glm::dvec3 intersect = glm::dvec3{0,0,0};
    glm::dvec3 normal = glm::dvec3{0,0,0};
    
    std::cout << "if intersection: " << s.intersect(r, intersect, normal) << std::endl;
    std::cout << "intersection point: " << glm::to_string(intersect) << std::endl;
    std::cout << "intersection normal: " << glm::to_string(normal) << std::endl;
};

void matrix_test() {
    glm::dvec3 a = glm::dvec3{1,0,1};
    glm::dvec3 b = glm::dvec3{0,2.5,0};
    glm::dvec3 c = glm::dvec3{3,3,3};
    std::cout << glm::to_string(glm::transpose(glm::mat3((a),(b),(c)))) << std::endl;
    std::cout << glm::dot(a,b) << std::endl;
//    Ray r = Ray(glm::dvec3{0,0,0}, glm::dvec3{0,2,2.5} - glm::dvec3{0,0,0});
//
//    ImpTriangle tri = ImpTriangle(glm::dvec3{1,2,2}, glm::dvec3{-1,2,2}, glm::dvec3{0,2,3});
//
//    glm::dvec3 intersect = glm::dvec3{0,0,0};
//    glm::dvec3 normal = glm::dvec3{0,0,0};
//
//    std::cout << "if intersection: " << tri.intersect(r, intersect, normal) << std::endl;
//    std::cout << "intersection point: " << glm::to_string(intersect) << std::endl;
//    std::cout << "intersection normal: " << glm::to_string(normal) << std::endl;
//
//    std::cout << "equal vec: " << glm::all(glm::equal(glm::dvec3{0,0,1}, glm::dvec3{0,0,0})) << std::endl;
};

void bbox_test() {
    BoundingBox b1 = BoundingBox(glm::dvec3{0,0,0}, glm::dvec3{2,2,2});
    BoundingBox b2 = BoundingBox(glm::dvec3{-2,-2,0}, glm::dvec3{1,1,2});
    glm::dvec3 point = glm::dvec3{1,1,1};
    
    std::cout << "if b1 and b2 intersect: " << b1.intersect(b2) << std::endl;
    std::cout << "point in b1: " << b1.contains(point) << std::endl;
};
