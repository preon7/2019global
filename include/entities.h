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

class ExpSphere : public Entity {
  
    int sectornum=10;
    int stacknum=10;
    // std::vector<float> vertices;
    float radius;
public:
    ExpSphere(glm::dvec3 pos, float radius) : Entity(),radius(radius) {
        this->pos = pos;
        std::vector<glm::vec3> vertices;
        std::vector<float> normals;
        std::vector<float> texCoords;

        float x, y, z,tmp;                              // vertex position
        // float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
        // float s, t;                                     // vertex texCoord

        float sectorStep = 2 * M_PI / sectornum;
        float stackStep = M_PI / stacknum;
        float sectorAngle, stackAngle;

        for(int i = 0; i <= stacknum; ++i)
        {
            stackAngle = M_PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
            tmp = radius * cosf(stackAngle);             // r * cos(u)
            z = radius * sinf(stackAngle) + pos.z;              // r * sin(u)

            // add (sectorCount+1) vertices per stack
            // the first and last vertices have same position and normal, but different tex coords
            for(int j = 0; j <= sectornum; ++j)
            {
                sectorAngle = j * sectorStep;           // starting from 0 to 2pi

                // vertex position (x, y, z)
                x = tmp * cosf(sectorAngle) + pos.x;             // r * cos(u) * cos(v)
                y = tmp * sinf(sectorAngle) + pos.y;
                
                vertices.push_back({x,y,z});

            }
        }
    }

    
  bool intersection(const Ray& ray,glm::dvec3 vertex0,glm::dvec3 vertex1,glm::dvec3 vertex2,
                             glm::dvec3& outIntersectionPoint)
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
          outIntersectionPoint = ray.origin + ray.dir * t;
          return true;
      }
      else // line intersection  not a ray intersection.
          return false;
  }
        
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(glm::vec3(this->pos.x - radius+0.1,this->pos.y - radius+0.1,this->pos.z - radius+0.1),
                                    glm::vec3(this->pos.x+0.1 + radius+0.1,this->pos.y + radius+0.1,this->pos.z + radius+0.1));
        return b;
    }

};


// TODO Implement explicit quad (triangles)
class ExpQuad : public Entity {
    std::vector<glm::vec3> vertices;
    //glm::vec3 vertices[];
    float width;
    float height;

public:
    ExpQuad(glm::dvec3 pos, float width, float height) : Entity(), width(width), height(height) {
        this->pos = pos;
        vertices.push_back({pos.x-width/2,pos.y-height/2,pos.z});

    }
    
    bool intersection(const Ray& ray,glm::dvec3 vertex0,glm::dvec3 vertex1,glm::dvec3 vertex2,
                               glm::dvec3& outIntersectionPoint)
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
            outIntersectionPoint = ray.origin + ray.dir * t;
            return true;
        }
        else // This means that there is a line intersection but not a ray intersection.
            return false;
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
    
    bool intersection(const Ray& ray,glm::dvec3 vertex0,glm::dvec3 vertex1,glm::dvec3 vertex2,
                               glm::dvec3& outIntersectionPoint)
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
            outIntersectionPoint = ray.origin + ray.dir * t;
            return true;
        }
        else // This means that there is a line intersection but not a ray intersection.
            return false;
    }
    
    
    
    BoundingBox boundingBox() const {
        BoundingBox b = BoundingBox(glm::vec3(pos.x- width/2,pos.y-length/2,pos.z-height/2),
                                    glm::vec3(pos.x+ width/2,pos.y+length/2,pos.z+height/2));
        return b;
    }

};
