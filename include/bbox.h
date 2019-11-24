#pragma once

#include <array>
#include <cassert>
#include <cstdlib>

#include <glm/glm.hpp>

/// Represents an axis-aligned bounding box.
struct BoundingBox {
    BoundingBox(glm::dvec3 min, glm::dvec3 max) : min(min), max(max) {
        assert(min.x < max.x);
        assert(min.y < max.y);
        assert(min.z < max.z);
    }

    double dx() const { return max.x - min.x; }
    double dy() const { return max.y - min.y; }
    double dz() const { return max.z - min.z; }

    const glm::dvec3 min;
    const glm::dvec3 max;

    /// Check if another bounding box intersects the current bounding box.
    bool intersect(const BoundingBox& other) const {
        // compute distance between two boxes
        glm::dvec3 pos_b1 = 0.5 * (min + max);
        glm::dvec3 pos_b2 = 0.5 * (other.min + other.max);
        
        glm::dvec3 distance = pos_b1 - pos_b2;
        
        // compare
        bool x_overlap = std::abs(distance.x) < (0.5*(max.x - min.x) + 0.5*(other.max.x - other.min.x));
        bool y_overlap = std::abs(distance.y) < (0.5*(max.y - min.y) + 0.5*(other.max.y - other.min.y));
        bool z_overlap = std::abs(distance.z) < (0.5*(max.z - min.z) + 0.5*(other.max.z - other.min.z));
        
        if (x_overlap && y_overlap && z_overlap) { return true; }
        else { return false; }
    }

    /// Check if a point lies within the bounding box.
    bool contains(glm::dvec3 point) const {
        glm::dvec3 pos = 0.5 * (min + max);
        
        bool x_in = std::abs(pos.x - point.x) <= 0.5*(max.x - min.x);
        bool y_in = std::abs(pos.y - point.y) <= 0.5*(max.y - min.y);
        bool z_in = std::abs(pos.z - point.z) <= 0.5*(max.z - min.z);
        
        if (x_in && y_in && z_in) { return true; }
        else { return false; }
    }
};
