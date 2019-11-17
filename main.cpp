#include <QApplication>

#include <iostream>

#include "camera.h"
#include "gui.h"

// for test
#include "ray.h"
#include "entities.h"
#include "glm/ext.hpp"

void test();

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    test();

    Camera camera({10, 0, 0});
    glm::dvec3 light{10, 10, 10};

    RayTracer raytracer(camera, light);

    // Set up scene
    Octree scene({-20, -20, -20}, {20, 20, 20});
    // TODO Add objects to the scene
    // scene.push_back(...);

    raytracer.setScene(&scene);

    Gui window(500, 500, raytracer);
    window.show();
    return app.exec();
}

void test() {
    ImpSphere s = ImpSphere(glm::dvec3{2,0,0}, 2);
    std::cout << 2 / 0 << std::endl;
    std::cout << "sphere created\n";
    std::cout << "sphere radius: " << s.radius << std::endl;
    std::cout << "sphere position: " << glm::to_string(s.pos) << std::endl;
    
    Ray r = Ray(glm::dvec3{0,0,0}, glm::dvec3{0,1,0});
    glm::dvec3 intersect = glm::dvec3{0,0,0};
    glm::dvec3 normal = glm::dvec3{0,0,0};
    
    std::cout << "if intersection: " << s.intersect(r, intersect, normal) << std::endl;
    std::cout << "intersection point: " << glm::to_string(intersect) << std::endl;
    std::cout << "intersection normal: " << glm::to_string(normal) << std::endl;
}
