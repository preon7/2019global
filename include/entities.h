#pragma once

#include <glm/glm.hpp>
#include <math.h>

#include "bbox.h"
#include "material.h"
#include "ray.h"

#include <iostream>
#include <algorithm>
#include <cassert>
#include <float.h>
#include "glm/ext.hpp"

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
    ImpSphere(glm::dvec3 pos, float radius, glm::dvec3 color) : Entity(Material(color)), radius(radius) {
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
            float v_1 = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
            float v_2 = (-b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
            
            // use the nearest one
            float base = std::min(std::abs(v_1), std::abs(v_2));
                intersect = glm::dvec3{base * a_1, base * a_2, base * a_3};
            
            // shift intersection point back
            intersect = intersect + ray.origin;
            normal = glm::normalize(intersect - this->pos);
            
            return true;
        };
        
        return false;
    }
    
    BoundingBox b = BoundingBox(glm::vec3(this->pos.x - radius,this->pos.y - radius,this->pos.z - radius),
                                glm::vec3(this->pos.x + radius,this->pos.y + radius,this->pos.z + radius));
    
    BoundingBox boundingBox() const {
        return b;
    }
};

// TODO Implement implicit triangle
class ImpTriangle : public Entity {
public:
    ImpTriangle(glm::dvec3 p1, glm::dvec3 p2, glm::dvec3 p3) : Entity(), p1(p1), p2(p2), p3(p3) {
        this->pos = 0.5*(0.5*(p1 + p2) + p3);
    }
    
    glm::dvec3 p1;
    glm::dvec3 p2;
    glm::dvec3 p3;
    
    glm::dvec3 edge1 = p2 - p1;
    glm::dvec3 edge2 = p3 - p1;
    glm::dvec3 normal = glm::normalize(glm::cross(edge1, edge2));
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        if (glm::dot(this->normal, ray.dir) == 0) { return false; }
        
        // point on plane: pos + (e1 * a + e2 * b) = ray.ori + c * ray.dir
        glm::mat3 A = glm::transpose(glm::mat3(edge1, edge2, -ray.dir));
        glm::mat3 A_i = glm::inverse(A);
        glm::dvec3 right = ray.origin - this->pos;
        glm::dvec3 sol = right * A_i;
        
//        std::cout << "pos: " << glm::to_string(this->pos) << std::endl;
//        std::cout << "edge1: " << glm::to_string(edge1) << std::endl;
//        std::cout << "edge2: " << glm::to_string(edge2) << std::endl;
//        std::cout << "inverse: " << glm::to_string(A_i) << std::endl;
//        std::cout << "right side: " << glm::to_string(right) << std::endl;
//        std::cout << "product: " << glm::to_string(sol) << std::endl;

        glm::dvec3 point = ray.origin + sol.z * ray.dir;
//        std::cout << "intermediate point: " << glm::to_string(point) << std::endl;
        
        if (glm::all(glm::equal(glm::normalize(point - p1), glm::normalize(p2 - point)))) {
            return true;
        }
        
        if (glm::all(glm::equal(glm::normalize(point - p2), glm::normalize(p3 - point)))) {
            return true;
        }
        
        if (glm::all(glm::equal(glm::normalize(point - p3), glm::normalize(p1 - point)))) {
            return true;
        }
        
        // position of point related to triangle
        glm::dvec3 d1 = glm::normalize(glm::cross(p1 - point, p2 - point));
        glm::dvec3 d2 = glm::normalize(glm::cross(p2 - point, p3 - point));
        glm::dvec3 d3 = glm::normalize(glm::cross(p3 - point, p1 - point));
        
        glm::dvec3 epsilon = glm::dvec3{1.0e-4, 1.0e-4, 1.0e-4};
        
//        std::cout << "d1: " << glm::to_string(d1) << std::endl;
//        std::cout << "d2: " << glm::to_string(d2) << std::endl;
//        std::cout << "d3: " << glm::to_string(d3) << std::endl;
//        std::cout << "eq1: " << glm::all(glm::lessThan(d1 - d2, epsilon)) << std::endl;
//        std::cout << "eq2: " << glm::all(glm::lessThan(d2 - d3, epsilon)) << std::endl;
        
        if (glm::all(glm::lessThan(d1 - d2, epsilon)) && glm::all(glm::lessThan(d2 - d3, epsilon))) {
            intersect = point;
            if (glm::dot(ray.dir, this->normal) < 0) {
                normal = this->normal;
            } else {
                normal = -this->normal;
            }
            
            return true;
        }
        
        return false;
    }
    
    glm::dvec3 min = glm::dvec3{std::min(std::min(p1.x, p2.x), p3.x),
        std::min(std::min(p1.y, p2.y), p3.y), std::min(std::min(p1.z, p2.z), p3.z)};
    glm::dvec3 max = glm::dvec3{std::max(std::max(p1.x, p2.x), p3.x),
        std::max(std::max(p1.y, p2.y), p3.y), std::max(std::max(p1.z, p2.z), p3.z)};
    
    BoundingBox b = BoundingBox(min, max);
    
    BoundingBox boundingBox() const {
        return b;
    }
};

class ExpRectangle : public Entity {
public:
    ExpRectangle(glm::dvec3 p1, glm::dvec3 p2, glm::dvec3 p3) : Entity(), p1(p1), p2(p2), p3(p3) {
        assert(glm::dot((p1-p3), (p2-p3)) == 0.0);
        this->pos = 0.5 * (p1 + p2);
    }
    
    glm::dvec3 p1;
    glm::dvec3 p2;  // p1, p2 is diagonal
    glm::dvec3 p3;
    
    glm::dvec3 p4 = this->pos + (this->pos - p3);  // p3, p4 is diagonal
    
    glm::dvec3 normal = glm::normalize(glm::cross((p1-p3), (p2-p3)));
    
    ImpTriangle t1 = ImpTriangle(p1, p2, p3);
    ImpTriangle t2 = ImpTriangle(p1, p2, p4);
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        if (t1.intersect(ray, intersect, normal)) {
            return true;
        }
        
        if (t2.intersect(ray, intersect, normal)) {
            return true;
        }
        
        return false;
    }
    
    glm::dvec3 min = glm::dvec3{std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)};
    glm::dvec3 max = glm::dvec3{std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)};
    BoundingBox b = BoundingBox(min, max);
    
    BoundingBox boundingBox() const {
        return b;
    }
};

class ExpBox : public Entity {
public:
    ExpBox(glm::dvec3 min, glm::dvec3 max) : Entity(), min(min), max(max) {
        assert(min.x < max.x);
        assert(min.y < max.y);
        assert(min.z < max.z);
    }
    
    const glm::dvec3 min;
    const glm::dvec3 max;
    glm::dvec3 down_left_bottom = min;
    glm::dvec3 down_right_bottom = glm::dvec3{max.x, min.y, min.z};
    glm::dvec3 down_left_top = glm::dvec3{min.x, max.y, min.z};
    glm::dvec3 down_right_top = glm::dvec3{max.x, max.y, min.z};
    
    glm::dvec3 up_left_bottom = glm::dvec3{min.x, min.y, max.z};
    glm::dvec3 up_right_bottom = glm::dvec3{max.x, min.y, max.z};
    glm::dvec3 up_left_top = glm::dvec3{min.x, max.y, max.z};
    glm::dvec3 up_right_top = max;
    
    std::array<std::unique_ptr<ExpRectangle>, 6> faces = {
        std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, up_right_bottom, up_left_bottom)),
        std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, up_left_top, down_left_top)),
        std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, down_right_top, down_left_top)),
        std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, up_left_bottom, up_left_top)),
        std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, down_right_bottom, down_right_top)),
        std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, down_left_top, down_right_top))
    };
    
//    faces[0] = std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, up_right_bottom, up_left_bottom));
//    faces[1] = std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, up_left_top, down_left_top));
//    faces[2] = std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, down_right_top, down_left_top));
//    faces[3] = std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, up_left_bottom, up_left_top));
//    faces[4] = std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, down_right_bottom, down_right_top));
//    faces[5] = std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, down_left_top, down_right_top));
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        
        
        //glm::dvec3 min_intersect = glm::dvec3{DBL_MAX, DBL_MAX, DBL_MAX};
        bool has_intersection = false;
        
        for (auto i = faces.begin(); i != faces.end(); ++i) {
            glm::dvec3 _intersect = {0,0,0};
            glm::dvec3 _normal = {0,0,0};
            double min_dist_square = DBL_MAX;
            
            if (i->get()->intersect(ray, _intersect, _normal)) {
                auto point_to = _intersect - ray.origin;
                double dist_square = pow(point_to.x, 2) + pow(point_to.y, 2) + pow(point_to.z, 2);
                
                if (dist_square < min_dist_square) {
                    min_dist_square = dist_square;
                    intersect = _intersect;
                    normal = _normal;
                }
                has_intersection = true;
            }
        }
            
        return has_intersection;
    }
    
    BoundingBox b = BoundingBox(min, max);
    
    BoundingBox boundingBox() const {
        return b;
    }
};

// TODO Implement explicit sphere (triangles)
// TODO Implement explicit quad (triangles)
// TODO Implement explicit cube (triangles)
