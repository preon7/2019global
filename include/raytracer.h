#pragma once

#include <algorithm>
#include <memory>
//#include <future>

#include <glm/glm.hpp>

#include "camera.h"
#include "entities.h"
#include "image.h"
#include "octree.h"
#include <float.h>

class RayTracer {
  public:
    RayTracer() = delete;
    RayTracer(const Camera& camera, glm::dvec3 light)
        : _camera(camera), _light(light), _image(std::make_shared<Image>(0,0)){};

    void setScene(const Octree* scene) { _scene = scene; }

    void run(int w, int h) {
        // TODO Implement this
        _image = std::make_shared<Image>(w, h);
        glm::dvec2 resolution = {0.0002, 0.0002};

        glm::dvec3 camrea_left = glm::normalize(glm::cross(_camera.up, _camera.forward));
        glm::dvec3 top_left = (_camera.pos + _camera.focalDist * _camera.forward +
                              camrea_left * double (w) * 0.5 * resolution.x + _camera.up * double (w) * 0.5 * resolution.y) - _camera.pos;
        // The structure of the for loop should remain for incremental rendering.
        for (int y = 0; y < h && _running; ++y) {
            for (int x = 0; x < w && _running; ++x) {
                // TODO Implement this
                // ray to the pixel
                glm::dvec3 direction = top_left - camrea_left * double (x) * resolution.x - _camera.up * double (y) * resolution.y;
//                std::cout << "ray to: " << glm::to_string(_camera.pos + direction) << std::endl;
                Ray r = Ray(_camera.pos, direction);
                
                std::vector<Entity*> objects = _scene->intersect(r);
                
                glm::dvec3 intersect = glm::dvec3{DBL_MAX, DBL_MAX, DBL_MAX};
                glm::dvec3 normal = glm::dvec3{0,0,0};
                
                Entity* front_obj;
                
                // check intersection
                for (int object_i = 0; object_i < objects.size(); object_i++) {
                    // std::cout << "checking intersection: " << object_i << std::endl;
                    glm::dvec3 current_intersect = glm::dvec3{0,0,0};
                    glm::dvec3 current_normal = glm::dvec3{0,0,0};
                    
                    // std::cout << "checking obj at: " << glm::to_string(objects[object_i]->pos) << std::endl;
                    if (objects[object_i]->intersect(r, current_intersect, current_normal)) {
                        if (glm::all(glm::lessThan(current_intersect - _camera.pos, intersect - _camera.pos))) {
                            intersect = current_intersect;
                            normal = current_normal;
                            front_obj = objects[object_i];
                            
//                            std::cout << "intersect at: " << glm::to_string(_camera.pos + r.dir) << std::endl;
                        }
                    }
                }
                
                if (front_obj) {
                    _image->setPixel(x, y, front_obj-> material.blinn_phong(r, _light, intersect, normal));
                } else {
                    _image->setPixel(x, y, {0, 0, 0});
                }
                
                front_obj = NULL;
            }
        }
    }

    bool running() const { return _running; }
    void stop() { _running = false; }
    void start() { _running = true; }

    std::shared_ptr<Image> getImage() const { return _image; }

  private:
    bool _running = false;
    const Octree* _scene;
    Camera _camera;
    glm::dvec3 _light;
    std::shared_ptr<Image> _image;
};
