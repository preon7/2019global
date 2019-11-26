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
void expsphere_test();
void boudingbox_test();

int main(int argc, char** argv) {
    
    QApplication app(argc, argv);
    entity_test();
    matrix_test();

    Camera camera({10, 0, 0});
    glm::dvec3 light{10, 10, 10};

    RayTracer raytracer(camera, light);
    

    // Set up scene
    Octree scene({-20, -20, -20}, {20, 20, 20});
    ImpSphere *s = new ImpSphere(glm::dvec3{2,0,0}, 2);
    
    scene.push_back(s);
    
    // TODO Add objects to the scene
    // scene.push_back(...);

    raytracer.setScene(&scene);
    
    expsphere_test();
    boudingbox_test();
    
    Gui window(500, 500, raytracer);
    
    window.show();
    // std::cout << "susscesssssssssssssss" << std::endl;
    return app.exec();
    

    

};

void entity_test() {
    ImpSphere s = ImpSphere(glm::dvec3{2,0,0}, 2);
    //std::cout << 2 / 0 << std::endl;
    std::cout << "sphere created" << std::endl;
    std::cout << "sphere radius: " << s.radius << std::endl;
    std::cout << "sphere position: " << glm::to_string(s.pos) << std::endl;
    
    Ray r = Ray(glm::dvec3{0,0,0}, glm::dvec3{0,1,0});
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


// jiaxin test

void expsphere_test() {
    ExpSphere expsphere = ExpSphere(glm::dvec3{1,0,0}, 10);
    //std::cout << 2 / 0 << std::endl;
    std::cout << "sphere created" << std::endl;
    std::cout << "sphere radius: " << expsphere.radius << std::endl;
    std::cout << "sphere position: " << glm::to_string(expsphere.pos) << std::endl;
    std::cout << "sphere position: " << glm::to_string(expsphere.pos) << std::endl;
    
    Ray r = Ray(glm::dvec3{1,2,3}, glm::dvec3{0,0,0});
    glm::dvec3 intersect = glm::dvec3{0,0,0};
    glm::dvec3 normal = glm::dvec3{0,0,0};
    
    std::cout << "if intersection: " << expsphere.intersect(r, intersect, normal) << std::endl;
    std::cout << "intersection point: " << glm::to_string(intersect) << std::endl;
    std::cout << "intersection normal: " << glm::to_string(normal) << std::endl;
}

void boudingbox_test() {
    BoundingBox b1 = BoundingBox(glm::dvec3{0,0,0}, glm::dvec3{2,2,2});
    BoundingBox b2 = BoundingBox(glm::dvec3{-2,-2,0}, glm::dvec3{1,1,2});
    glm::dvec3 point = glm::dvec3{1,1,1};
    
    std::cout << "if b1 and b2 intersect: " << b1.intersect(b2) << std::endl;
    std::cout << "point in b1: " << b1.contains(point) << std::endl;
};
