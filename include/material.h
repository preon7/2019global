#pragma once

#include <glm/glm.hpp>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "ray.h"

/// Represents the material properties of an entity. For now it only contains color, but it should
/// probably be extended to allow more options.
struct Material {
    constexpr explicit Material(glm::dvec3 color) : color(std::move(color)) {
        specular_color = {1,1,1};
        diffuse_color = color * 0.5;
    }
    
    Material(glm::dvec3 color, glm::dvec3 shader) : color(std::move(color)), shader_parameters(shader) {
        specular_color = {1,1,1};
        diffuse_color = color * 0.5;
    }

    glm::dvec3 color;
    glm::dvec3 diffuse_color;
    glm::dvec3 specular_color;
    
    glm::dvec3 shader_parameters = glm::dvec3(0.1,0.7,1);
    
    glm::dvec3 blinn_phong(Ray ray, glm::dvec3 light, glm::dvec3 intersect, glm::dvec3 normal) {
        double specular_power = 5;
        
        glm::dvec3 la = color * shader_parameters.x;
        glm::dvec3 ld = std::max(0.0, glm::dot(normal, glm::normalize(light - intersect))) * diffuse_color * shader_parameters.y;
        glm::dvec3 bisector = glm::normalize(glm::normalize(-ray.dir) + glm::normalize(light - intersect));
        glm::dvec3 ls = pow(std::max(0.0, glm::dot(normal, bisector)), specular_power) * specular_color * shader_parameters.z;
       
        auto output = la + ld + ls;
        
//        double max_val = std::max(std::max(output.x, output.y), output.z);
//        if (max_val > 1.0) { output = output / max_val; }
        
        // prevent the value from exceeding 1
        output = {std::min(output.x, 1.0), std::min(output.y, 1.0), std::min(output.z, 1.0)};
        
        return output;
    }
};
