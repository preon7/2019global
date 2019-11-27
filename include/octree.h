#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "bbox.h"
#include "entities.h"

class Octree {
  public:
    Octree(glm::dvec3 min, glm::dvec3 max) : min(min), max(max), _root(Node({min, max})) {}
    
    glm::dvec3 min;
    glm::dvec3 max;

    /// Store an entity in the correct position of the octree.
    void push_back(Entity* object) {
        // TODO Implement this
        if (!_root._bbox.intersect(object->boundingBox())) {
            return;  // raise error
        }
        _root.push_obj(object);
        //Node *n = _root._children[0].get();
        
        // _root.partition();
//        for(auto i = _root._children.begin(); i != _root._children.end(); ++i) {
//            if (i->get()->_bbox.intersect(object->boundingBox())) {
//                i->get()->_entities.push_back(object);
//            }
//        }
    }

    /// Returns list of entities that have the possibility to be intersected by the ray.
    std::vector<Entity*> intersect(const Ray& ray) const {
        // TODO Implement this
        //return _root._entities;
//        std::vector<Entity*> out_list;
//
//        for (auto i = _root._children.begin(); i != _root._children.end(); ++i) {
//            ExpBox intersect_box = ExpBox(i->get()->_bbox.min, i->get()->_bbox.max);
//
//            glm::dvec3 intersect = glm::dvec3{0,0,0};
//            glm::dvec3 normal = glm::dvec3{0,0,0};
//
//            if (intersect_box.intersect(ray, intersect, normal)) {
//                // append entities in child to output list
//                // TODO: remove dup
//                out_list.insert(out_list.end(), i->get()->_entities.begin(), i->get()->_entities.end());
//            }
//        }
//
        return _root.intersect(ray);
    }

  private:
    struct Node {
        explicit Node(const BoundingBox& bbox) : _bbox(bbox) {}

        /// Subdivides the current node into 8 children.
        void partition(){
            if (!is_leaf()) {
                return;
            }
            
            glm::dvec3 middle = (_bbox.min + _bbox.max) * 0.5;
            
            // check if all 8 children intersect with all objects
            bool all_in = true;
            for (int e_i = 0; e_i < _entities.size(); e_i++) {
                BoundingBox _b = _entities[e_i]->boundingBox();
                bool b_l_corner = glm::all(glm::lessThan(_bbox.min, _b.min)) && glm::all(glm::lessThan(_b.min, middle));
                bool t_r_corner = glm::all(glm::lessThan(middle, _b.max)) && glm::all(glm::lessThan(_b.max, _bbox.max));
                all_in = all_in && b_l_corner && t_r_corner;
            }
            
            if (all_in) { return; }
            
            // split into eight equal boxes
            _children[0] = std::make_unique<Node>(Node(BoundingBox(_bbox.min, middle)));
            _children[1] = std::make_unique<Node>(Node(BoundingBox(glm::dvec3{_bbox.min.x, middle.y, _bbox.min.z},
                                                                   glm::dvec3{middle.x, _bbox.max.y, middle.z})));
            _children[2] = std::make_unique<Node>(Node(BoundingBox(glm::dvec3{middle.x, _bbox.min.y, _bbox.min.z},
                                                                   glm::dvec3{_bbox.max.x, middle.y, middle.z})));
            _children[3] = std::make_unique<Node>(Node(BoundingBox(glm::dvec3{middle.x, middle.y, _bbox.min.z},
                                                                   glm::dvec3{_bbox.max.x, _bbox.max.y , middle.z})));
            
            _children[4] = std::make_unique<Node>(Node(BoundingBox(middle, _bbox.max)));
            _children[5] = std::make_unique<Node>(Node(BoundingBox(glm::dvec3{_bbox.min.x, middle.y, middle.z},
                                                                   glm::dvec3{middle.x, _bbox.max.y, _bbox.max.z})));
            _children[6] = std::make_unique<Node>(Node(BoundingBox(glm::dvec3{middle.x, _bbox.min.y, middle.z},
                                                                   glm::dvec3{_bbox.max.x, middle.y, _bbox.max.z})));
            _children[7] = std::make_unique<Node>(Node(BoundingBox(glm::dvec3{_bbox.min.x, _bbox.min.y, middle.z},
                                                                   glm::dvec3{middle.x, middle.y , _bbox.max.z})));
            
        };

        bool is_leaf() const { return _children[0] == nullptr; }
        
        // propagate object deep into tree
        void push_obj(Entity* object) {
            _entities.push_back(object);
            partition();
            if (is_leaf()) {
                return;
            }
            
            for (int c_i = 0; c_i < 8; c_i++) {
                if (glm::all(glm::lessThan(_children[c_i]->_bbox.min, object->boundingBox().min)) && glm::all(glm::lessThan(object->boundingBox().max, _children[c_i]->_bbox.max))) {
                    _children[c_i]->push_obj(object);
                }
            }
        }
        
        std::vector<Entity*> intersect(Ray ray) const {
            if (is_leaf()) {
                return _entities;
            }
            
            std::vector<Entity*> out_list;
            
            for (auto i = _children.begin(); i != _children.end(); ++i) {
                ExpBox intersect_box = ExpBox(i->get()->_bbox.min, i->get()->_bbox.max);
                
                glm::dvec3 intersect = glm::dvec3{0,0,0};
                glm::dvec3 normal = glm::dvec3{0,0,0};
                
                if (intersect_box.intersect(ray, intersect, normal)) {
                    // append entities in child to output list
                    // TODO: remove dup
                    out_list.insert(out_list.end(), i->get()->intersect(ray).begin(), i->get()->intersect(ray).end());
                }
            }
            
            return out_list;
        }

        BoundingBox _bbox;
        std::vector<Entity*> _entities;
        std::array<std::unique_ptr<Node>, 8> _children;
    };

    Node _root;
};
