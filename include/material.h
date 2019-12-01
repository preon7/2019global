#pragma once

#include <glm/glm.hpp>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "ray.h"

/// Represents the material properties of an entity. For now it only contains color, but it should
/// probably be extended to allow more options.
struct Material {
    constexpr explicit Material(glm::dvec3 color) : color(std::move(color)) {}

    glm::dvec3 color;
    glm::dvec3 shader_parameters = glm::dvec3(0.05,0.3,1);
    
    glm::dvec3 blinn_phong(Ray ray, glm::dvec3 light, glm::dvec3 intersect, glm::dvec3 normal) {
        double specular_power = 5;
        
        glm::dvec3 specular_color = {1,1,1};
        
        glm::dvec3 la = color * shader_parameters.x;
        glm::dvec3 ld = std::max(0.0, glm::dot(normal, glm::normalize(light - intersect))) * color * shader_parameters.y;
        glm::dvec3 bisector = glm::normalize(glm::normalize(-ray.dir) + glm::normalize(light - intersect));
        glm::dvec3 ls = pow(std::max(0.0, glm::dot(normal, bisector)), specular_power) * specular_color * shader_parameters.z;
        
        auto output = la +ld +ls;
        output = {std::min(output.x, 1.0), std::min(output.y, 1.0), std::min(output.z, 1.0)};
        
        return output;
    }
};
