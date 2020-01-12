#pragma once

#include <algorithm>
#include <memory>
//#include <future>

#include <glm/glm.hpp>

#include "camera.h"
#include "entities.h"
#include "image.h"
#include "octree.h"
#include "ray.h"
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
                
                if (y == 3 && x == w-1) {
                    int stop = 1;
                }
                
                glm::dvec3 direction = top_left - camrea_left * double (x) * resolution.x - _camera.up * double (y) * resolution.y;
                //std::cout << "ray to: " << glm::to_string(_camera.pos + direction) << std::endl;
                Ray r = Ray(_camera.pos, direction);
                
                std::vector<Entity*> objects = _scene->intersect(r);
                
                glm::dvec3 intersect = glm::dvec3{DBL_MAX, DBL_MAX, DBL_MAX};
                glm::dvec3 normal = glm::dvec3{0,0,0};
                
//                Entity* front_obj = nullptr;
//
//                // check intersection
//                for (int object_i = 0; object_i < objects.size(); object_i++) {
//                    // std::cout << "checking intersection: " << object_i << std::endl;
//                    glm::dvec3 current_intersect = glm::dvec3{0,0,0};
//                    glm::dvec3 current_normal = glm::dvec3{0,0,0};
//
//                    double min_dist_square = DBL_MAX;
//
//                    // std::cout << "checking obj at: " << glm::to_string(objects[object_i]->pos) << std::endl;
//                    if (objects[object_i]->intersect(r, current_intersect, current_normal)) {
//                        auto point_to = current_intersect - r.origin;
//                        double dist_square = pow(point_to.x, 2) + pow(point_to.y, 2) + pow(point_to.z, 2);
//
//                        if (dist_square < min_dist_square) {
//                            min_dist_square = dist_square;
//                            intersect = current_intersect;
//                            normal = current_normal;
//                            front_obj = objects[object_i];
//
//                            //std::cout << "intersect at: " << glm::to_string(_camera.pos + r.dir) << std::endl;
//                        }
//                    }
//                }
                Entity* front_obj = get_front_obj(objects, r, intersect, normal);
                
                if (front_obj) {
                    auto coord = front_obj->getTextureCoord(intersect);
                    _image->setPixel(x, y, front_obj->material.blinn_phong_texture(r, _light, intersect, normal, std::get<0>(coord), std::get<1>(coord)));
                    //_image->setPixel(x, y, front_obj->material.blinn_phong(r, _light, intersect, normal));
                    
//                    // check visibility
//                    glm::dvec3 float_intersect = {0,0,0};
//                    if (glm::dot(r.dir, normal) < 0)
//                        float_intersect = intersect + _epsilon * normal;
//                    else
//                        float_intersect = intersect - _epsilon * normal;
//
//                    Ray r_light = Ray(float_intersect, _light - float_intersect);
//                    auto block_obj = get_front_obj(objects, r_light, intersect, normal);

//                    // set color
//                    if (block_obj)
//                        _image->setPixel(x, y, {0,0,0});
//                    else
//                        _image->setPixel(x, y, reflect_color(objects, front_obj, r, intersect, normal, 1));
                        _image->setPixel(x, y, refraction_color(objects, front_obj, r, intersect, normal, 1));
                } else {
                    _image->setPixel(x, y, {0, 0, 0});
                }
                
                front_obj = NULL;
            }
        }
    }
    
    Entity* get_front_obj(std::vector<Entity*> objects, Ray r, glm::dvec3& intersect, glm::dvec3& normal) {
        Entity* front_obj = nullptr;
        
        // check intersection
        for (int object_i = 0; object_i < objects.size(); object_i++) {
            // std::cout << "checking intersection: " << object_i << std::endl;
            glm::dvec3 current_intersect = glm::dvec3{0,0,0};
            glm::dvec3 current_normal = glm::dvec3{0,0,0};
            
            double min_dist_square = DBL_MAX;
            
            // std::cout << "checking obj at: " << glm::to_string(objects[object_i]->pos) << std::endl;
            if (objects[object_i]->intersect(r, current_intersect, current_normal)) {
                auto point_to = current_intersect - r.origin;
                double dist_square = pow(point_to.x, 2) + pow(point_to.y, 2) + pow(point_to.z, 2);
                
                if (dist_square < min_dist_square) {
                    min_dist_square = dist_square;
                    intersect = current_intersect;
                    normal = current_normal;
                    front_obj = objects[object_i];
                    
                    //std::cout << "intersect at: " << glm::to_string(_camera.pos + r.dir) << std::endl;
                }
            }
        }
        
        return front_obj;
    }
    
    int max_reflect = 2;
    double reflection_rate = 0.3;
    
    glm::dvec3 reflect_color(std::vector<Entity*> objects, Entity* front_obj, Ray r, glm::dvec3 intersect, glm::dvec3 normal, int reflect_time) {
        if (reflect_time > max_reflect) {
            return {0,0,0};
        }
        
        if (glm::dot(r.dir, normal) > 0) {
            normal = -normal;
        }
        
        auto coord = front_obj->getTextureCoord(intersect);
        auto local_color = front_obj->material.blinn_phong_texture(r, _light, intersect, normal, std::get<0>(coord), std::get<1>(coord));
        
        //double epsilon = 1e-3;
        auto float_intersect = intersect + _epsilon * normal;  // chose a point float on the surface to avoid self-collision
        auto reflected_ray = Ray(float_intersect,
                                 2. * normal + r.dir); // 2*normal - (-dir)
        
        auto new_front = get_front_obj(objects, reflected_ray, float_intersect, normal);
        
        glm::dvec3 reflected_color;
        if (new_front) {
            reflected_color = reflect_color(objects, new_front, reflected_ray, intersect, normal, reflect_time + 1);
        } else {
            reflected_color = glm::dvec3{0,0,0};
        }
        
        auto output = local_color + reflection_rate * reflected_color;
        
        return {std::min(output.x, 1.0), std::min(output.y, 1.0), std::min(output.z, 1.0)};
    }
    
 
    glm::dvec3 refraction_color(std::vector<Entity*> objects, Entity* front_obj, Ray r, glm::dvec3 intersect, glm::dvec3 normal, int times) {
        
        if (times > 2) {
            return {0,0,0};
        }

        double eta_i = 1.0;
        double eta_t = front_obj->getRefractIndex();
//        material.refraction;
        double transparancy = front_obj->getTransparancy();
        double kr = fresnel(r.dir, normal, eta_t);
        std::cout << kr <<std::endl;
        double i_dot_n = glm::dot(r.dir, normal);
        if (i_dot_n < 0.0) {
            //Outside the surface
            i_dot_n = -i_dot_n;
        } else {
            //Inside the surface; invert the normal and swap the indices of refraction
            normal = -normal;
            eta_i = eta_t;
            eta_t = 1.0;
        }
        double eta = eta_i / eta_t;
        double k = 1.0 - (eta * eta) * (1.0 - i_dot_n * i_dot_n);
        
//        double kr = fresnel(r.dir, normal, eta_t);
//        double transparancy = 0.5;
   
        auto coord = front_obj->getTextureCoord(intersect);
        auto local_color = front_obj->material.blinn_phong_texture(r, _light, intersect, normal, std::get<0>(coord), std::get<1>(coord));
        
        auto float_intersect = intersect + _epsilon * normal;  // chose a point float on the surface to avoid self-collision
        auto reflected_ray = Ray(float_intersect, 2. * normal + r.dir);
//                                 r.dir - 2.0*glm::dot(normal, r.dir)* normal );
        auto refract_dir = (r.dir + i_dot_n *normal)*eta - normal * sqrt(k);
        auto refracted_ray = Ray(float_intersect, refract_dir);
        
        auto new_front1 = get_front_obj(objects, refracted_ray, float_intersect, normal);
        auto new_front2 = get_front_obj(objects, reflected_ray, float_intersect, normal);
        
        glm::dvec3 reflected_color;
        glm::dvec3 refracted_color;
        if (new_front1) {
            refracted_color = refraction_color(objects,new_front1, refracted_ray, intersect, normal, times+1);
        } else {
            refracted_color = glm::dvec3{0,0,0};
        }

        if (new_front2) {
            reflected_color = reflect_color(objects,new_front2, reflected_ray, intersect, normal, times+1);
        } else {
            refracted_color = glm::dvec3{0,0,0};
        }
        
        auto output = reflected_color*kr + refracted_color*(1-kr) ;
        output = local_color * transparancy + output;//* transparancy;
//        local_color + reflected_color *(1 - transparancy) + refracted_color * transparancy;
               
        return {std::min(output.x, 1.0), std::min(output.y, 1.0), std::min(output.z, 1.0)};
    }
    
    
    double fresnel(glm::dvec3 incident,glm::dvec3 normal,double index)  {
        double i_dot_n =  glm::dot(incident,normal);
        double eta_i = 1.0;
        double eta_t = index;
        if (i_dot_n > 0.0) {
            eta_i = eta_t;
            eta_t = 1.0;
        }

        double sin_t = eta_i / eta_t * sqrtf(std::max( 0.0,(1.0 - i_dot_n * i_dot_n)));
        if(sin_t >= 1.0) {
            //Total internal reflection
            return 1.0;
        } else {
            double cos_t = sqrtf (std::max(0.0,(1.0 - sin_t * sin_t)));
            double cos_i = abs(i_dot_n);
            double r_s = ((eta_t * cos_i) - (eta_i * cos_t)) / ((eta_t * cos_i) + (eta_i * cos_t));
            double r_p = ((eta_i * cos_i) - (eta_t * cos_t)) / ((eta_i * cos_i) + (eta_t * cos_t));
            return (r_s * r_s + r_p * r_p) / 2.0;
        }
    }
    

    bool running() const { return _running; }
    void stop() { _running = false; }
    void start() { _running = true; }

    std::shared_ptr<Image> getImage() const { return _image; }

  private:
    double _epsilon = 1e-3;
    bool _running = false;
    const Octree* _scene;
    Camera _camera;
    glm::dvec3 _light;
    std::shared_ptr<Image> _image;
};
