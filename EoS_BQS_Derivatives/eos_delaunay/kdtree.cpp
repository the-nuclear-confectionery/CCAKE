#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
 
#include "kdtree.h"
 
template<typename coordinate_type, size_t dimensions>
std::ostream& operator<<(std::ostream& out, const point<coordinate_type, dimensions>& pt) {
    out << '(';
    for (size_t i = 0; i < dimensions; ++i) {
        if (i > 0)
            out << ", ";
        out << pt.get(i);
    }
    out << ')';
    return out;
}

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
	template<typename coordinate_type, size_t dimensions>
	kdtree<coordinate_type, dimensions>::node*
	kdtree<coordinate_type, dimensions>::make_tree(size_t begin, size_t end, size_t index) {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        auto i = nodes_.begin();
        std::nth_element(i + begin, i + n, i + end, node_cmp(index));
        index = (index + 1) % dimensions;
        nodes_[n].left_ = make_tree(begin, n, index);
        nodes_[n].right_ = make_tree(n + 1, end, index);
        return &nodes_[n];
    }
 
	template<typename coordinate_type, size_t dimensions>
    void kdtree::nearest(kdtree::node* root, const point_type& point, size_t index) {
        if (root == nullptr)
            return;
        ++visited_;
        double d = root->distance(point);
        if (best_ == nullptr || d < best_dist_) {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index);
    }

    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    kdtree::kdtree(iterator begin, iterator end) : nodes_(begin, end) {
		for ( size_t inode = 0; inode < nodes_.size(); inode++ )
			nodes_[inode].oindex_ = inode; // added by me
        root_ = make_tree(0, nodes_.size(), 0);
    }
 
    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool kdtree::empty() const { return nodes_.empty(); }
 
    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t kdtree::visited() const { return visited_; }
 
    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double kdtree::distance() const { return std::sqrt(best_dist_); }
 
    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& kdtree::nearest(const point_type& pt ) {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& kdtree::nearest(const point_type& pt, size_t & best_oindex ) {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
		best_oindex = best_->oindex_;
        return best_->point_;
    }

