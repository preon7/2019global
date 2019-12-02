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

#define PI 3.1415926535

/// A base class for all entities in the scene.
struct Entity {

    constexpr Entity() : material(Material(glm::dvec3(1, 0, 0))) {}
    constexpr Entity(const Material& material) : material(material) {}

    /// Check if a ray intersects the object. The arguments intersect and normal will contain the
    /// point of intersection and its normals.
    virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const = 0;

    /// Returns an axis-aligned bounding box of the entity.
    virtual BoundingBox boundingBox() const = 0;
    
    // returns the texture coordinates
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const = 0;

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
    
    glm::dvec3 up_vec = glm::dvec3{0,0,radius};
    //glm::dvec3 left_vec = pos + glm::dvec3{0,radius,0};
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        double unit_length_v = 2.0 * PI * radius / 320.0; // 10 repeating patterns on equator
        glm::dvec3 to_inter = intersect - pos;
        
        //double cos_hori = glm::dot(glm::dvec3{intersect.x, intersect.y, pos.z}, left_vec) / pow(radius, 2);
        double cos_vert = glm::dot(to_inter, up_vec) / pow(radius, 2);
        double angle_to_up = acos(cos_vert);
        
        // calculate relative y coordinate
        int y = int ((0.5 * PI * radius - radius * angle_to_up) / unit_length_v);
        
        double small_r = radius * sin(angle_to_up);
        glm::dvec3 left_middle_vec = glm::dvec3{0,small_r,0};
        
        double cos_hori = glm::dot(glm::dvec3{to_inter.x, to_inter.y, 0}, left_middle_vec) / pow(small_r, 2);
        
        double unit_length_h = 2.0 * PI * small_r / 320.0; // 10 repeating patterns
        
        // calculate relative x coordinate
        int x = int (small_r * acos(cos_hori) / unit_length_h);
        
        return std::make_tuple(x,y);
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
        
//        if (glm::all(glm::equal(glm::normalize(point - p1), glm::normalize(p2 - point)))) {
//            return true;
//        }
//
//        if (glm::all(glm::equal(glm::normalize(point - p2), glm::normalize(p3 - point)))) {
//            return true;
//        }
//
//        if (glm::all(glm::equal(glm::normalize(point - p3), glm::normalize(p1 - point)))) {
//            return true;
//        }
        
        // position of point related to triangle
        glm::dvec3 d1 = glm::normalize(glm::cross(p1 - point, p2 - point));
        glm::dvec3 d2 = glm::normalize(glm::cross(p2 - point, p3 - point));
        glm::dvec3 d3 = glm::normalize(glm::cross(p3 - point, p1 - point));
        
//        std::cout << "1-p: " << glm::to_string(p1 - point) << std::endl;
//        std::cout << "2-p: " << glm::to_string(p2 - point) << std::endl;
//        std::cout << "3-p: " << glm::to_string(p3 - point) << std::endl;
        
//        glm::dvec3 epsilon = glm::dvec3{1.0e-3, 1.0e-3, 1.0e-3};
        
//        std::cout << "d1: " << glm::to_string(d1) << std::endl;
//        std::cout << "d2: " << glm::to_string(d2) << std::endl;
//        std::cout << "d3: " << glm::to_string(d3) << std::endl;
//        std::cout << "eq1: " << glm::all(glm::lessThan(d1 - d2, epsilon)) << std::endl;
//        std::cout << "eq2: " << glm::all(glm::lessThan(d2 - d3, epsilon)) << std::endl;
        double epsilon = 1.0e-3;
        auto diff_1 = d1 - d2;
        auto diff_2 = d2 - d3;
        bool cp_1 = pow(diff_1.x, 2) + pow(diff_1.y, 2) + pow(diff_1.z, 2) < epsilon;
        bool cp_2 = pow(diff_2.x, 2) + pow(diff_2.y, 2) + pow(diff_2.z, 2) < epsilon;
        
        if (cp_1 && cp_2) {  // glm::all(glm::lessThan(d1 - d2, epsilon)) && glm::all(glm::lessThan(d2 - d3, epsilon))
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
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
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
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
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
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
    }
};

// TODO Implement explicit sphere (triangles)
// TODO Implement explicit quad (triangles)
// TODO Implement explicit cube (triangles)

class ExpSphere : public Entity {
    
public:
    
    ExpSphere(glm::dvec3 pos, float radius, glm::dvec3 color) : Entity(Material(color)), radius(radius){
        this->pos = pos;
        float x, y, z,tmp;     // vertex position
        float sectorStep = 2 * M_PI / sectornum;
        float stackStep = M_PI / stacknum;
        float sectorAngle, stackAngle;
        // stack angle = pi/2 - pi*stackstep/stacknum
        //sector angle = 2*pi *sectorstep/sectornum
        for(int i = 0; i <= stacknum; ++i)  // interate the surface vertically
        {
            stackAngle = M_PI / 2 - i * stackStep;  // starting from pi/2 to -pi/2
            
            tmp = radius * cos(stackAngle);    // r * cos(sectorAngle)
            z = radius * sin(stackAngle);  // z = r * sin(sectorAngle)
            
            for(int j = 0; j <= sectornum; ++j) // interate the surface horizontally
            {
                sectorAngle = j * sectorStep;
                //std::cout << "sectorAngle is:" << sectorAngle << std::endl;
                x = tmp * cos(sectorAngle);    // x = r * cos(sectorAngle)*cos(sectorAngle)
                y = tmp * sin(sectorAngle);    // y = r*cos(sectorAngle)*sin(sectorAngle)
                glm::dvec3 vertex = glm::dvec3{x,y,z} - pos;
                vertices.push_back(vertex);
                normal.push_back(glm::normalize(vertex));
                // std::cout << "normal vertices is:" << glm::to_string(normal) << std::endl;
                
            }
        }
        
        // triangulate adjacent vertices to form polygons
        int k1, k2; // two adjencent vertices
        for(int i = 0; i < stacknum; ++i) // interate all trangles vertically
        {
            k1 = i * (sectornum + 1);     // current stack
            k2 = k1 + sectornum + 1;      // next stack
            for(int j = 0; j < sectornum; ++j, ++k1, ++k2)
            {
                if(i != 0)
                    triangles.push_back(new ImpTriangle(vertices.at(k1),vertices.at(k2),vertices.at(k1+1)));
                if(i != (stacknum-1))
                    triangles.push_back(new ImpTriangle(vertices.at(k1+1),vertices.at(k2),vertices.at(k2+1)));
                //                std::cout<<"trangle push back success"<<std::endl;
            }
        }
    }
    float radius;
    int sectornum=10; // number of sectors on the surface (horizontal)
    int stacknum=10; // number of stacks on the surface (vertical)
    std::vector<glm::dvec3> vertices;
    std::vector<glm::dvec3> normal;
    std::vector<Entity*> triangles;
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        bool flag = false;
        
        double min_dist_sqare = DBL_MAX;
        glm::dvec3 min_intersect = glm::dvec3{DBL_MAX, DBL_MAX, DBL_MAX};
        glm::dvec3 current_normal = {0,0,0};
        for (int i = 1; i < triangles.size(); i++) {
            //            std::cout << "trangle size is:" << triangles.size() << std::endl;
            if (triangles[i]->intersect(ray, intersect, normal) == true) {
                glm::dvec3 to_point = intersect - ray.origin;
                if (pow(to_point.x,2) + pow(to_point.y,2) + pow(to_point.z,2) <= min_dist_sqare ) {
                    min_intersect = intersect;
                    current_normal = normal;
                    min_dist_sqare = pow(to_point.x,2) + pow(to_point.y,2) + pow(to_point.z,2);
                }
                flag = true;
            }
        }
        
        normal = current_normal;
        intersect = min_intersect;
        return flag;
    }
    
    
    BoundingBox b = BoundingBox(glm::vec3(this->pos.x - radius-1,this->pos.y - radius-1,this->pos.z - radius-1),
                                glm::vec3(this->pos.x + radius+1,this->pos.y + radius+1,this->pos.z + radius+1));
    
    BoundingBox boundingBox() const {
        return b;
    }
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
    }
};


// TODO Implement explicit quad (triangles)
class ExpQuad : public Entity {
public:
    ExpQuad(glm::dvec3 pos, float width, float length) : Entity(), width(width), length(length) {
        this->pos = pos;
        vertices.push_back(glm::dvec3{pos.x+width/2,pos.y+length/2,pos.z}); //upright
        vertices.push_back(glm::dvec3{pos.x-width/2,pos.y+length/2,pos.z}); //upleft
        vertices.push_back(glm::dvec3{pos.x+width/2,pos.y-length/2,pos.z}); //downright
        vertices.push_back(glm::dvec3{pos.x-width/2,pos.y-length/2,pos.z}); //downleft
        //        std::cout << "point: " << glm::to_string(vertices.at(0)) << std::endl;
        //        std::cout << "point: " << glm::to_string(vertices.at(1)) << std::endl;
        //        std::cout << "point: " << glm::to_string(vertices.at(2)) << std::endl;
        //        std::cout << "point: " << glm::to_string(vertices.at(3)) << std::endl;
        
    }
    std::vector<glm::vec3> vertices;
    float width;
    float length;
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        // test if the ray intersects with any trangle on the surface
        bool flag = false;
        std::cout<<"ray.org is:"<<glm::to_string(ray.origin)<<std::endl;
        std::cout<<"ray.dir is:"<<glm::to_string(ray.dir)<<std::endl;
        
        if(ray_trangle_intersection(ray, intersect,vertices.at(1), vertices.at(2), vertices.at(0)))
            flag = true;
        if(ray_trangle_intersection(ray, intersect,vertices.at(1), vertices.at(3), vertices.at(2)))
            flag = true;
        
        if(flag == true)
            return true;
        else
            return false;
    }
    
    bool ray_trangle_intersection(const Ray& ray, glm::dvec3& intersect, glm::dvec3 vertex0, glm::dvec3 vertex1, glm::dvec3 vertex2) const{
        
        const float EPSILON = 0.0000001f; //threshold
        glm::dvec3 edge1, edge2, N;
        // compute plane's normal
        
        edge1 = vertex1 - vertex0;
        edge2 = vertex2 - vertex0;
        N = glm::cross(edge1,edge2); // normal vector that perpendicular to the plane
        std::cout << "N1: " << glm::to_string(N) << std::endl;
        //step 1: find intersection point p=O+tR; t distance from ray origin O to p
        //test if ray and plane are parallel
        float N_RayDir = glm::dot(N, ray.dir);
        std::cout << "N_RayDir: " << N_RayDir << std::endl;
        if(fabs(N_RayDir)< EPSILON)
            return false;
        //if the ray direction is perpendicular(dotproduct=0) to the N, ray is parallel to the plane
        
        float d = glm::dot(N, vertex0); // Ax+By+Cz+D=0 D is the distance from the origin to the plane
        std::cout << "d: " << d << std::endl;
        //compute t: distance from ray origin O to p
        float t = (d - glm::dot(N, ray.origin) )/N_RayDir;
        std::cout << "t: " << t << std::endl;
        if(t < 0) return false; //trangle is behind the ray
        // compute the intersection point
        intersect = ray.origin + double (t) * ray.dir;
        std::cout << "intersect!: " << glm::to_string(intersect) << std::endl;
        
        // step 2: if P is inside the triangle or not
        glm::dvec3 c;
        // test edge 0
        glm::dvec3 e0 = vertex1 - vertex0;
        glm::dvec3 pv0 = intersect - vertex0;
        c = glm::cross(e0, pv0);
        if(glm::dot(N,c) < 0) return false;
        
        // test edge 1
        glm::dvec3 e1 = vertex2 - vertex1;
        glm::dvec3 pv1 = intersect - vertex1;
        c = glm::cross(e1, pv1);
        if(glm::dot(N,c) < 0) return false;
        
        // test edge 2
        glm::dvec3 e2 = vertex0 - vertex2;
        glm::dvec3 pv2 = intersect - vertex2;
        c = glm::cross(e2, pv2);
        if(glm::dot(N,c) < 0) return false;
        
        return true; // ray hits inside the trangle
    }
    
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(glm::vec3(pos.x- width/2,pos.y-length/2,pos.z),
                                    glm::vec3(pos.x+ width/2,pos.y+length/2,pos.z));
        return b;
    }
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
    }
};

// TODO Implement explicit cube (triangles)
class ExpCube : public Entity {
    std::vector<glm::vec3> vertices;
    float width, length, height;
public:
    ExpCube(glm::dvec3 pos, float width, float length, float height) : Entity(), width(width), length(length), height(height)
    {
        this->pos = pos;
        vertices.push_back({pos.x-width/2,pos.y-length/2,pos.z-height/2});
        vertices.push_back({pos.x-width/2,pos.y-length/2,pos.z+height/2});
        vertices.push_back({pos.x+width/2,pos.y-length/2,pos.z-height/2});
        vertices.push_back({pos.x+width/2,pos.y-length/2,pos.z+height/2});
        vertices.push_back({pos.x-width/2,pos.y+length/2,pos.z+height/2});
        vertices.push_back({pos.x-width/2,pos.y+length/2,pos.z-height/2});
        vertices.push_back({pos.x+width/2,pos.y+length/2,pos.z-height/2});
        vertices.push_back({pos.x+width/2,pos.y+length/2,pos.z+height/2});
        
        
    }
    
    //    bool ray_trangle_intersection(const Ray& ray,glm::dvec3 vertex0,glm::dvec3 vertex1,glm::dvec3 vertex2, glm::dvec3& intersect) const
    //    {
    //        const float EPSILON = 0.0000001;
    //        glm::dvec3 edge1, edge2, h, s, q;
    //        float a,f,u,v;
    //        edge1 = vertex1 - vertex0;
    //        edge2 = vertex2 - vertex0;
    //        h = glm::cross(ray.dir,edge2);
    //        a = glm::dot(edge1,h);
    //        if (a > -EPSILON && a < EPSILON)
    //            return false;    // This ray is parallel to this triangle.
    //        f = 1.0/a;
    //        s = ray.origin - vertex0;
    //        u = f * glm::dot(s,h);
    //        if (u < 0.0 || u > 1.0)
    //            return false;
    //        q = glm::cross(s,edge1);
    //        v = f * glm::dot(ray.dir,q);
    //        if (v < 0.0 || u + v > 1.0)
    //            return false;
    //        // At this stage we can compute t to find out where the intersection point is on the line.
    //        double t = f * glm::dot(edge2,q);
    //        if (t > EPSILON && t < 1/EPSILON) // ray intersection
    //        {
    //            intersect = ray.origin + ray.dir * t;
    //            return true;
    //        }
    //        else // This means that there is a line intersection but not a ray intersection.
    //            return false;
    //    }
    
    bool intersection(const Ray& ray, std::vector<glm::vec3> vertices) const {
        
        // test if the ray intersects with any trangle on the surface
        glm::dvec3 intersect;
        bool flag = false;
        if(ray_trangle_intersection(ray,intersect,vertices.at(0),vertices.at(1),vertices.at(2)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(3),vertices.at(1),vertices.at(2)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(4),vertices.at(5),vertices.at(7)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(7),vertices.at(5),vertices.at(6)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(1),vertices.at(0),vertices.at(4)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(4),vertices.at(0),vertices.at(5)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(3),vertices.at(7),vertices.at(2)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(7),vertices.at(6),vertices.at(2)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(1),vertices.at(4),vertices.at(3)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(3),vertices.at(4),vertices.at(7)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(0),vertices.at(5),vertices.at(2)))
            flag = true;
        if(ray_trangle_intersection(ray,intersect,vertices.at(2),vertices.at(5),vertices.at(6)))
            flag = true;
        
        if(flag == true)
            return true;
        else
            return false;
    }
    
    bool ray_trangle_intersection(const Ray& ray, glm::dvec3& intersect, glm::dvec3 vertex0, glm::dvec3 vertex1, glm::dvec3 vertex2) const{
        
        const float EPSILON = 0.0000001f; //threshold
        glm::dvec3 edge1, edge2, N;
        // compute plane's normal
        
        edge1 = vertex1 - vertex0;
        edge2 = vertex2 - vertex0;
        N = glm::cross(edge1,edge2); // vector that perpendicular to the plane
        
        //step 1: find intersection point p=O+tR; t distance from ray origin O to p
        //test if ray and plane are parallel
        float N_RayDir = glm::dot(N, ray.dir);
        if(fabs(N_RayDir)< EPSILON)
            return false;
        //if the ray direction is perpendicular(dotproduct=0) to the N, ray is parallel to the plane
        
        float d = glm::dot(N, vertex0); // Ax+By+Cz+D=0 D is the distance from the origin to the plane
        
        //compute t: distance from ray origin O to p
        float t = -(glm::dot(N, ray.origin) + d)/N_RayDir;
        if(t < 0) return false; //trangle is behind the ray
        // compute the intersection point
        intersect = ray.origin + double (t) * ray.dir;
        // std::cout << "t: " << t << std::endl;
        
        // step 2: if P is inside the triangle or not
        glm::dvec3 c;
        // test edge 0
        glm::dvec3 e0 = vertex1 - vertex0;
        glm::dvec3 pv0 = intersect - vertex0;
        c = glm::cross(e0, pv0);
        if(glm::dot(N,c) < 0) return false;
        
        // test edge 1
        glm::dvec3 e1 = vertex2 - vertex1;
        glm::dvec3 pv1 = intersect - vertex1;
        c = glm::cross(e1, pv1);
        if(glm::dot(N,c) < 0) return false;
        
        // test edge 2
        glm::dvec3 e2 = vertex0 - vertex2;
        glm::dvec3 pv2 = intersect - vertex2;
        c = glm::cross(e2, pv2);
        if(glm::dot(N,c) < 0) return false;
        
        return true; // ray hits inside the trangle
    }
    
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(glm::vec3(pos.x- width/2,pos.y-length/2,pos.z-height/2),
                                    glm::vec3(pos.x+ width/2,pos.y+length/2,pos.z+height/2));
        return b;
    }
    
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
    }
};

// TODO Implement explicit cone (triangles)
class ExpCone : public Entity {
public:
    ExpCone(glm::dvec3 pos, float height, float radius, glm::dvec3 color): Entity(Material(color)), height(height), radius(radius){
        this->pos = pos; // Top vertex of cone
        float alpha;
        glm::dvec3 loc = glm::dvec3{0,0,0};
        vertices.push_back(pos);
        double numSubdivisions = 50.0;
        for (int i = 0; i <= numSubdivisions; ++i) {
            alpha = i * 360.0/numSubdivisions;
            loc.x = pos.x + radius * cos(alpha * PI / 180.0);
            //            std::cout << "alpha: " << alpha << std::endl;
            //            std::cout << "cos(alpha): " << cos(alpha * PI / 180.0) << std::endl;
            loc.y = pos.y + radius * sin(alpha * PI / 180.0);
            loc.z = pos.z - height;
            vertices.push_back(loc);
        }
        
        for(int i = 1; i < vertices.size()-1; ++i){
            triangles.push_back(new ImpTriangle(vertices.at(0),vertices.at(i),vertices.at(i+1)));
        }
    }
    std::vector<glm::dvec3> vertices; // all vertices of cone
    std::vector<Entity*> triangles;
    float height, radius;
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        // test if the ray intersects with any trangle on the surface
        bool flag = false;
        //triangles.push_back(new ImpTriangle(vertices.at(0),vertices.at(numtri-1),vertices.at(1)));
        
        double min_dist_sqare = DBL_MAX;
        glm::dvec3 min_intersect = glm::dvec3{DBL_MAX, DBL_MAX, DBL_MAX};
        glm::dvec3 current_normal = {0,0,0};
        for (int i = 0; i < triangles.size(); i++) {
            if (triangles[i]->intersect(ray, intersect, normal) == true) {
                glm::dvec3 to_point = intersect - ray.origin;
                
                if (pow(to_point.x,2) + pow(to_point.y,2) + pow(to_point.z,2) <= min_dist_sqare ) {
                    min_intersect = intersect;
                    current_normal = normal;
                    min_dist_sqare = pow(to_point.x,2) + pow(to_point.y,2) + pow(to_point.z,2);
                }
                
                flag = true;
            }
        }
        
        normal = current_normal;
        intersect = min_intersect;
        
        return flag;
    }
    
    
    
    BoundingBox b = BoundingBox(glm::vec3(this->pos.x- radius,this->pos.y - radius,this->pos.z-height),
                                glm::vec3(this->pos.x+ radius,this->pos.y + radius,this->pos.z));
    BoundingBox boundingBox() const {
        
        return b;
    }
  
    virtual std::tuple<int,int> getTextureCoord(glm::dvec3 intersect) const {
        //TODO: implement
        return std::make_tuple(0,0);
    }
};
