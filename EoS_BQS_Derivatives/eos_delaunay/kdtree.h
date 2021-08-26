#ifndef KDTREE_H
#define KDTREE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
 
template<typename coordinate_type, size_t dimensions>
class point {
public:
    point(std::array<coordinate_type, dimensions> c) : coords_(c) {}
    point(std::initializer_list<coordinate_type> list);
    coordinate_type get(size_t index) const;
    double distance(const point& pt) const;
private:
    std::array<coordinate_type, dimensions> coords_;
};
 
template<typename coordinate_type, size_t dimensions>
std::ostream& operator<<(std::ostream& out, const point<coordinate_type, dimensions>& pt);
 
/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
template<typename coordinate_type, size_t dimensions>
class kdtree {
public:
    typedef point<coordinate_type, dimensions> point_type;
private:
    struct node {
        node(const point_type& pt) : point_(pt), left_(nullptr), right_(nullptr) {}
        coordinate_type get(size_t index) const {
            return point_.get(index);
        }
        double distance(const point_type& pt) const {
            return point_.distance(pt);
        }
        point_type point_;
        node* left_;
        node* right_;
		size_t oindex_;	// added by me
    };
    node* root_ = nullptr;
    node* best_ = nullptr;
    double best_dist_ = 0;
    size_t visited_ = 0;
    std::vector<node> nodes_;
 
    struct node_cmp {
        node_cmp(size_t index) : index_(index) {}
        bool operator()(const node& n1, const node& n2) const {
            return n1.point_.get(index_) < n2.point_.get(index_);
        }
        size_t index_;
    };
 
    node* make_tree(size_t begin, size_t end, size_t index);
 
    void nearest(node* root, const point_type& point, size_t index);
public:
    kdtree(const kdtree&) = delete;
    kdtree& operator=(const kdtree&) = delete;
 
    template<typename iterator>
    kdtree(iterator begin, iterator end);

	kdtree(){};	// add a default constructor (added by me)
 
    //template<typename func>
    //kdtree(func&& f, size_t n);
 
    bool empty() const;
    size_t visited() const;
    double distance() const;
    const point_type& nearest(const point_type& pt );
    const point_type& nearest(const point_type& pt, size_t & best_oindex );
};


#endif 
