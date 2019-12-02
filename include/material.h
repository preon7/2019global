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
    
    double specular_power = 5;
    
    glm::dvec3 blinn_phong(Ray ray, glm::dvec3 light, glm::dvec3 intersect, glm::dvec3 normal) {
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
    
    glm::dvec3 blinn_phong_texture(Ray ray, glm::dvec3 light, glm::dvec3 intersect, glm::dvec3 normal, int relative_x, int relative_y) {
        Texture t = Texture(color);
        glm::dvec3 texture_color = t.color_at(relative_x, relative_y);
        glm::dvec3 texture_diffuse_color = texture_color * 0.5;
        
        glm::dvec3 la = texture_color * shader_parameters.x;
        glm::dvec3 ld = std::max(0.0, glm::dot(normal, glm::normalize(light - intersect))) * texture_diffuse_color * shader_parameters.y;
        glm::dvec3 bisector = glm::normalize(glm::normalize(-ray.dir) + glm::normalize(light - intersect));
        glm::dvec3 ls = pow(std::max(0.0, glm::dot(normal, bisector)), specular_power) * specular_color * shader_parameters.z;
        
        auto output = la + ld + ls;
        output = {std::min(output.x, 1.0), std::min(output.y, 1.0), std::min(output.z, 1.0)};
        
        return output;
    }
    
private:
    struct Texture {
        explicit Texture(glm::dvec3 color) {
//            this->width = 32;
//            this->height = 32;
            
            //int pattern[width][height][3];
            for (int x=0; x<width; x++) {
                for (int y=0; y<height; y++) {
                    if (x <= int (width / 2) && y <= int (height / 2)) {
                        pattern[x][y][0] = 1;
                        pattern[x][y][1] = 1;
                        pattern[x][y][2] = 1;
                        continue;
                    }
                    
                    if (x > int (width / 2) && y > int (height / 2)) {
                        pattern[x][y][0] = 1;
                        pattern[x][y][1] = 1;
                        pattern[x][y][2] = 1;
                        continue;
                    }
                    
                    pattern[x][y][0] = color.x;
                    pattern[x][y][1] = color.y;
                    pattern[x][y][2] = color.z;
                }
            }
        }
        
        static const int width = 32;
        static const int height = 32;
        int pattern[width][height][3];
        
        glm::dvec3 color_at(int x, int y) {
            glm::dvec3 output = {pattern[x % width][y % height][0],pattern[x%width][y%height][1],pattern[x%width][y%height][2]};
            
            return output;
        }
    };
};
