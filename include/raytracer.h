#pragma once

#include <algorithm>
#include <memory>
//#include <future>

#include <glm/glm.hpp>

#include "camera.h"
#include "entities.h"
#include "image.h"
#include "octree.h"

class RayTracer {
  public:
    RayTracer() = delete;
    RayTracer(const Camera& camera, glm::dvec3 light)
        : _camera(camera), _light(light), _image(std::make_shared<Image>(0,0)){};

    void setScene(const Octree* scene) { _scene = scene; }

    void run(int w, int h) {
        // TODO Implement this
        _image = std::make_shared<Image>(w, h);
        
        glm::dvec3 u=glm::cross(_camera.up, _camera.forward);
        u = glm::normalize(u);
        glm::dvec3 v =glm::cross(_camera.forward, u);
        v = glm::normalize(v);
        
        glm::dvec3 c = _camera.pos- _camera.forward * _camera.focalDist;
        // c = glm::normalize(c);
        glm::dvec3 l = c - multiply(u, w/2)- multiply(v, h/2);
        // l = glm::normalize(l);
        
        // The structure of the for loop should remain for incremental rendering.
        for (int y = 0; y < h && _running; ++y) {
            for (int x = 0; x < w && _running; ++x) {
                // TODO Implement this
                //  _image->setPixel(x, y, {0, 0, 0});
                
                double factor_x = double (x) * double (w) / 32.0;
                double factor_y = double (y) * double (h) / 32.0;
                glm::dvec3 imgloc = l + u * factor_x + v * factor_y;
                // imgloc = glm::normalize(imgloc);
                std::cout << "img location: " << glm::to_string(imgloc) << std::endl;
    
                Ray ray = Ray(_camera.pos, imgloc);
                std::vector<Entity*> objlist = _scene->intersect(ray);
                float minDistance = 1000000;
                for(int i=0; i< objlist.size(); ++i){
                    float dist = distance(objlist.at(i),_camera.pos);
                    if(dist < minDistance){
                        minDistance = dist;
                        auto obj = objlist.at(i)->pos;
                        _image->setPixel(x,y,{obj.r,obj.g,obj.b});
//                        std::cout <<"  x="<<x<<"  y="<<y<< "  i="<<i<< std::endl;
//                        std::cout <<"  r="<<obj.r<<"  g="<<obj.g<< "  b="<<obj.b<< std::endl;
//                        std::cout<<" dist2 ="<<dist<<std::endl;
                    }
                }

                
    
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
    
    glm::dvec3 multiply(glm::dvec3 vec, double n){
        glm::dvec3 new_vec;
        new_vec.x=vec.x*n;
        new_vec.y=vec.y*n;
        new_vec.z=vec.z*n;
        return new_vec;
    };
    
    glm::dvec3 vecadd(glm::dvec3 vec1, glm::vec3 vec2){
        glm::dvec3 new_vec;
        new_vec.x=vec1.x + vec2.x;
        new_vec.y=vec1.y + vec2.y;
        new_vec.z=vec1.z + vec2.z;
        return new_vec;
    };
    
    float distance(Entity* v1, glm::dvec3 v2){
        return sqrt(pow(v1->pos.x-v2.x,2)+pow(v1->pos.y-v2.y,2)+pow(v1->pos.z-v2.z,2));
    }
};




