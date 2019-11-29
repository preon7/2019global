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
class ImpTriangle : Entity {
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
        
        glm::dvec3 epsilon = glm::dvec3{1.0e-5, 1.0e-5, 1.0e-5};
        
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
    
    BoundingBox boundingBox() const {
        glm::dvec3 min = glm::dvec3{std::min(std::min(p1.x, p2.x), p3.x),
            std::min(std::min(p1.y, p2.y), p3.y), std::min(std::min(p1.z, p2.z), p3.z)};
        glm::dvec3 max = glm::dvec3{std::max(std::max(p1.x, p2.x), p3.x),
            std::max(std::max(p1.y, p2.y), p3.y), std::max(std::max(p1.z, p2.z), p3.z)};
        
        BoundingBox b = BoundingBox(min, max);
        return b;
    }
};

class ExpRectangle : Entity {
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
    
    BoundingBox boundingBox() const {
        glm::dvec3 min = glm::dvec3{std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)};
        glm::dvec3 max = glm::dvec3{std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)};
        BoundingBox b = BoundingBox(min, max);
        return b;
    }
};

class ExpBox : Entity {
public:
    ExpBox(glm::dvec3 min, glm::dvec3 max) : Entity(), min(min), max(max) {
        assert(min.x < max.x);
        assert(min.y < max.y);
        assert(min.z < max.z);
    }
    
    const glm::dvec3 min;
    const glm::dvec3 max;
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        glm::dvec3 down_left_bottom = min;
        glm::dvec3 down_right_bottom = glm::dvec3{max.x, min.y, min.z};
        glm::dvec3 down_left_top = glm::dvec3{min.x, max.y, min.z};
        glm::dvec3 down_right_top = glm::dvec3{max.x, max.y, min.z};
        
        glm::dvec3 up_left_bottom = glm::dvec3{min.x, min.y, max.z};
        glm::dvec3 up_right_bottom = glm::dvec3{max.x, min.y, max.z};
        glm::dvec3 up_left_top = glm::dvec3{min.x, max.y, max.z};
        glm::dvec3 up_right_top = max;
        
        std::array<std::unique_ptr<ExpRectangle>, 6> faces;
        
        faces[0] = std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, up_right_bottom, up_left_bottom));
        faces[1] = std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, up_left_top, down_left_top));
        faces[2] = std::make_unique<ExpRectangle>(ExpRectangle(down_left_bottom, down_right_top, down_left_top));
        faces[3] = std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, up_left_bottom, up_left_top));
        faces[4] = std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, down_right_bottom, down_right_top));
        faces[5] = std::make_unique<ExpRectangle>(ExpRectangle(up_right_top, down_left_top, down_right_top));
        
        glm::dvec3 min_intersect = glm::dvec3{DBL_MAX, DBL_MAX, DBL_MAX};
        bool has_intersection = false;
        
        for (auto i = faces.begin(); i != faces.end(); ++i) {
            glm::dvec3 _intersect = {0,0,0};
            glm::dvec3 _normal = {0,0,0};
            if (i->get()->intersect(ray, _intersect, _normal)) {
                if (glm::all(glm::lessThan(_intersect, min_intersect))) {
                    min_intersect = _intersect;
                    normal = _normal;
                    intersect = _intersect;
                }
                has_intersection = true;
            }
        }
            
        return has_intersection;
    }
    
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(min, max);
        return b;
    }
};

// TODO Implement explicit sphere (triangles)

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
                //std::cout << "vertices is:" << glm::to_string(vertex) << std::endl;
                vertices.push_back(vertex);
                normal.push_back(glm::normalize(vertex));
                // std::cout << "normal vertices is:" << glm::to_string(normal) << std::endl;
                
            }
        }
    }
    float radius;
    int sectornum=10; // number of sectors on the surface (horizontal)
    int stacknum=10; // number of stacks on the surface (vertical)
    std::vector<glm::dvec3> vertices;
    std::vector<glm::dvec3> normal;
    
    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        std::vector<int> indices;
        // triangulate adjacent vertices to form polygons
        int k1, k2; // two adjencent vertices
        for(int i = 0; i < stacknum; ++i) // interate all trangles vertically
        {
            k1 = i * (sectornum + 1);     // current stack
            k2 = k1 + sectornum + 1;      // next stack
            // k1 - k1+1
            //  |  /
            //  k2

            for(int j = 0; j < sectornum; ++j, ++k1, ++k2)
            {
                if(i != 0)
                {
                    indices.push_back(k1);
                    indices.push_back(k2);
                    indices.push_back(k1 + 1);
                }
                if(i != (stacknum-1))
                {
                    indices.push_back(k1 + 1);
                    indices.push_back(k2);
                    indices.push_back(k2 + 1);
                }
            }
        }
        // test if the ray intersects with any trangle on the surface
        
        bool flag = false;
        
        for(int i = 0; i< indices.size(); i=i+3){
            
            bool label = ray_trangle_intersection(ray,intersect, vertices.at(indices.at(i)), vertices.at(indices.at(i+1)), vertices.at(indices.at(i+2)));
            // std::cout << "point: " << glm::to_string(intersect) << std::endl;
            if(label== true)
                flag = true;

        }
        normal = glm::normalize(intersect);
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
        BoundingBox b = BoundingBox(glm::vec3(this->pos.x - radius+0.1,this->pos.y - radius+0.1,this->pos.z - radius+0.1),
                                    glm::vec3(this->pos.x+0.1 + radius+0.1,this->pos.y + radius+0.1,this->pos.z + radius+0.1));
        return b;
    }

};


// TODO Implement explicit quad (triangles)
class ExpQuad : public Entity {
public:
    ExpQuad(glm::dvec3 pos, float width, float height) : Entity(), width(width), height(height) {
        this->pos = pos;
        vertices.push_back(glm::dvec3{pos.x+width/2,pos.y+height/2,pos.z}); //upright
        vertices.push_back(glm::dvec3{pos.x-width/2,pos.y+height/2,pos.z}); //upleft
        vertices.push_back(glm::dvec3{pos.x+width/2,pos.y-height/2,pos.z}); //downright
        vertices.push_back(glm::dvec3{pos.x-width/2,pos.y-height/2,pos.z}); //downleft
//        std::cout << "point: " << glm::to_string(vertices.at(0)) << std::endl;
//        std::cout << "point: " << glm::to_string(vertices.at(1)) << std::endl;
//        std::cout << "point: " << glm::to_string(vertices.at(2)) << std::endl;
//        std::cout << "point: " << glm::to_string(vertices.at(3)) << std::endl;

    }
    std::vector<glm::vec3> vertices;
    float width;
    float height;

    bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
        // test if the ray intersects with any trangle on the surface
        bool flag = false;
        std::cout<<"ray.org is:"<<glm::to_string(ray.origin)<<std::endl;
        std::cout<<"ray.dir is:"<<glm::to_string(ray.dir)<<std::endl;
        
        if(ray_trangle_intersection(ray, intersect,vertices.at(0), vertices.at(1), vertices.at(2)))
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
          
          //step 1: find intersection point p=O+tR; t distance from ray origin O to p
          //test if ray and plane are parallel
          float N_RayDir = glm::dot(N, ray.dir);
          if(fabs(N_RayDir)< EPSILON)
              return false;
          //if the ray direction is perpendicular(dotproduct=0) to the N, ray is parallel to the plane
          
          float d = glm::dot(N, vertex1); // Ax+By+Cz+D=0 D is the distance from the origin to the plane
          
          //compute t: distance from ray origin O to p
          float t = (glm::dot(N, ray.origin) + d)/N_RayDir;
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
        BoundingBox b = BoundingBox(glm::vec3(pos.x- width/2,pos.y-height/2,pos.z),
                                    glm::vec3(pos.x+ width/2,pos.y+height/2,pos.z));
        return b;
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
    
    bool ray_trangle_intersection(const Ray& ray,glm::dvec3 vertex0,glm::dvec3 vertex1,glm::dvec3 vertex2, glm::dvec3& intersect) const
    {
        const float EPSILON = 0.0000001;
        glm::dvec3 edge1, edge2, h, s, q;
        float a,f,u,v;
        edge1 = vertex1 - vertex0;
        edge2 = vertex2 - vertex0;
        h = glm::cross(ray.dir,edge2);
        a = glm::dot(edge1,h);
        if (a > -EPSILON && a < EPSILON)
            return false;    // This ray is parallel to this triangle.
        f = 1.0/a;
        s = ray.origin - vertex0;
        u = f * glm::dot(s,h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = glm::cross(s,edge1);
        v = f * glm::dot(ray.dir,q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        double t = f * glm::dot(edge2,q);
        if (t > EPSILON && t < 1/EPSILON) // ray intersection
        {
            intersect = ray.origin + ray.dir * t;
            return true;
        }
        else // This means that there is a line intersection but not a ray intersection.
            return false;
    }
    
    bool intersection(const Ray& ray, std::vector<glm::vec3> vertices) const {

        // test if the ray intersects with any trangle on the surface
        glm::dvec3 intersect;
        bool flag = false;
        if(ray_trangle_intersection(ray,vertices[0],vertices[1],vertices[2],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[3],vertices[1],vertices[2],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[4],vertices[5],vertices[7],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[7],vertices[5],vertices[6],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[1],vertices[0],vertices[4],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[4],vertices[0],vertices[5],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[3],vertices[7],vertices[2],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[7],vertices[6],vertices[2],intersect))
            flag = true;
        if(ray_trangle_intersection(ray,vertices[1],vertices[4],vertices[3],intersect))
             flag = true;
        if(ray_trangle_intersection(ray,vertices[3],vertices[4],vertices[7],intersect))
             flag = true;
        if(ray_trangle_intersection(ray,vertices[0],vertices[5],vertices[2],intersect))
             flag = true;
        if(ray_trangle_intersection(ray,vertices[2],vertices[5],vertices[6],intersect))
             flag = true;
        
        if(flag == true)
            return true;
        else
            return false;
    }
    
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(glm::vec3(pos.x- width/2,pos.y-length/2,pos.z-height/2),
                                    glm::vec3(pos.x+ width/2,pos.y+length/2,pos.z+height/2));
        return b;
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
            int numSubdivisions = 10;
            for (int i = 0; i < numSubdivisions; ++i) {
                alpha = i * 360/numSubdivisions;
                loc.x = pos.x + radius * cos(alpha);
                loc.y = pos.y + radius * sin(alpha);
                loc.z = pos.z - height;
                vertices.push_back(loc);
            }
        }
        std::vector<glm::dvec3> vertices; // all vertices of cone
        float height, radius;
     
        bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const {
                   // test if the ray intersects with any trangle on the surface
            bool flag = false;
            int numtri = vertices.size();
            if(ray_trangle_intersection(ray,intersect,vertices.at(0),vertices.at(numtri-1),vertices.at(1))) flag = true;
            for(int i = 0; i < numtri-2; ++i){
                if(ray_trangle_intersection(ray,intersect,vertices.at(0),vertices.at(i),vertices.at(i+1))){
                    std::cout << "intersection!!!: " << glm::to_string(intersect) << std::endl;
                    flag = true;
                }
                    
            }
            
            normal = glm::normalize(intersect);
        
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
              float t = (glm::dot(N, ray.origin) + d)/N_RayDir;
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
            BoundingBox b = BoundingBox(glm::vec3(pos.x- radius,pos.y-radius,pos.z), glm::vec3(pos.x+ radius,pos.y+radius,pos.z-height));
            return b;
        }
    
};
