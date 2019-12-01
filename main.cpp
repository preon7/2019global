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
void Exp_test();

int main(int argc, char** argv) {
    QApplication app(argc, argv);
//    entity_test();

    //Camera camera({10, 0, 0});
    Camera camera({-10, 0, 0}, {1, 0, 0}, 0.1);
    glm::dvec3 light{-10, 10, 10};

    RayTracer raytracer(camera, light);

//    Exp_test();

    // Set up scene
    Octree scene({-20, -20, -20}, {20, 20, 20});
//    ImpSphere *s = new ImpSphere(glm::dvec3{4,0,0}, 2, {0,1,0});
    ImpSphere *s2 = new ImpSphere(glm::dvec3{3,4,4}, 2, {1,0,0});
    ImpSphere *s3 = new ImpSphere(glm::dvec3{4,-4,4}, 2, {0,0,1});
    ExpCone *cone = new ExpCone( glm::dvec3{0,0,2}, 4 ,2, {1,0,1} );
//    ImpTriangle *t = new ImpTriangle({1,0,3}, {0,-1,0}, {0,3,0});
//    ExpSphere *es = new ExpSphere( glm::dvec3{4,2,2}, 2, {1,0,1} ); // 2,2,1


//    scene.push_back(s);
    scene.push_back(s2);
    scene.push_back(s3);
    scene.push_back(cone);
//    scene.push_back(es);
    // TODO Add objects to the scene
    // scene.push_back(...);

    raytracer.setScene(&scene);

    Gui window(500, 500, raytracer);
    window.show();
    return app.exec();
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

//jiaxin test
void Exp_test() {
    //sphere
    ExpSphere s = ExpSphere(glm::dvec3{0,0,0}, 5 ,glm::dvec3{1,1,1});
    std::cout << "ExpSphere created" << std::endl;
    std::cout << "ExpSphere radius: " << s.radius << std::endl;
    std::cout << "ExpSphere position: " << glm::to_string(s.pos) << std::endl;
    
    //quad
//    ExpQuad s = ExpQuad(glm::dvec3{1,1,1}, 10 ,10);
//    std::cout << "ExpQuad created" << std::endl;
//    std::cout << "ExpQuad width: " << s.width << std::endl;
//    std::cout << "ExpQuad position: " << glm::to_string(s.pos) << std::endl;
    
    //cone
//    ExpCone s = ExpCone( glm::dvec3{1,1,1}, 10 ,10, glm::dvec3{1,0,0} );
//    std::cout << "ExpCone created" << std::endl;
//    std::cout << "ExpCone radius: " << s.radius << std::endl;
//    std::cout << "ExpCone position: " << glm::to_string(s.pos) << std::endl;
    
    Ray r = Ray(glm::dvec3{10,0,0}, glm::dvec3{1,0,0});
    std::cout << "ray direction: " << glm::to_string(r.dir) << std::endl;
    glm::dvec3 intersect = glm::dvec3{0,0,0};
    glm::dvec3 normal = glm::dvec3{0,0,0};
    
    std::cout << "if intersection: " << s.intersect(r, intersect, normal) << std::endl;
    std::cout << "intersection point: " << glm::to_string(intersect) << std::endl;
    std::cout << "intersection normal: " << glm::to_string(normal) << std::endl;
};
