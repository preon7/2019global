#pragma once

#include <glm/glm.hpp>
#include <math.h>

#include "bbox.h"
#include "material.h"
#include "ray.h"

#include <iostream>
#include <algorithm>

/// A base class for all entities in the scene.
struct Entity {

    constexpr Entity() : material(Material(glm::dvec3(1, 0, 0))) {}
    constexpr Entity(const Material& material) : material(material) {}

    /// Check if a ray intersects the object. The arguments intersect and normal will contain the
    /// point of intersection and its normals.
    virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const = 0;

    /// Returns an axis-aligned bounding box of the entity.
    virtual BoundingBox boundingBox() const = 0;

    glm::dvec3 pos = {0, 0, 0};
    Material material;
};

// TODO Implement implicit sphere
class ImpSphere : public Entity {
public:
    ImpSphere(glm::dvec3 pos, float radius) : Entity(), radius(radius) {
        this->pos = pos;
    }
    
    float radius;
    
    // (x-p)^2 + (y-p)^2 + (z-p)^2 = r^2
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        // shift the position
        glm::dvec3 new_pos = this->pos - ray.origin;
        
        float a_1 = 1;
        float a_2 = 1;
        float a_3 = 1;
        if (ray.dir.x != 0) {
            a_2 = ray.dir.y / ray.dir.x;
            a_3 = ray.dir.z / ray.dir.x;
        } else if (ray.dir.y != 0) {
            a_1 = ray.dir.x / ray.dir.y;
            a_3 = ray.dir.z / ray.dir.y;
        } else if (ray.dir.z != 0) {
            a_2 = ray.dir.y / ray.dir.z;
            a_1 = ray.dir.x / ray.dir.z;
        } else { return false; }
        
        float a = pow(a_1, 2) + pow(a_2, 2) + pow(a_3, 2);
        float b = -2 * (new_pos.x * a_1 + new_pos.y * a_2 + new_pos.z * a_3);
        float c = pow(new_pos.x, 2) + pow(new_pos.y, 2) + pow(new_pos.z, 2) - pow(radius, 2);
        
        if (pow(b, 2) - 4 * a * c < 0) {
//            std::cout << a << std::endl;
//            std::cout << b << std::endl;
//            std::cout << c << std::endl;
            return false;
        } else {
            float v_1 = (-b + sqrt(pow(b, 2) - 4 * a * c)) / 2 * a;
            float v_2 = (-b - sqrt(pow(b, 2) - 4 * a * c)) / 2 * a;
            
            // use the nearest one
            float base = std::min(v_1, v_2);
                intersect = glm::dvec3{base * a_1, base * a_2, base * a_3};
            
            // shift intersection point back
            intersect = intersect + ray.origin;
            normal = glm::normalize(intersect - this->pos);
            
            return true;
        };
        
        return false;
    }
    
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(glm::vec3(this->pos.x - radius,this->pos.y - radius,this->pos.z - radius),
                                    glm::vec3(this->pos.x + radius,this->pos.y + radius,this->pos.z + radius));
        return b;
    }
};
// TODO Implement implicit triangle

// TODO Implement explicit sphere (triangles)
// TODO Implement explicit quad (triangles)
// TODO Implement explicit cube (triangles)
