#ifndef KDTREE_H
#define KDTREE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
 
/**
 * Class for representing a point. coordinate_type must be a numeric type.
 */
template<typename coordinate_type, size_t dimensions>
class point {
public:
    point(std::array<coordinate_type, dimensions> c) : coords_(c) {}
    point(std::initializer_list<coordinate_type> list) {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }
    /**
     * Returns the coordinate in the given dimension.
     *
     * @param index dimension index (zero based)
     * @return coordinate in the given dimension
     */
    coordinate_type get(size_t index) const {
        return coords_[index];
    }
    /**
     * Returns the distance squared from this point to another
     * point.
     *
     * @param pt another point
     * @return distance squared from this point to the other point
     */
    double distance(const point& pt) const {
        double dist = 0;
        for (size_t i = 0; i < dimensions; ++i) {
            double d = get(i) - pt.get(i);
            dist += d * d;
        }
        return dist;
    }
    double logdistance(const point& pt) const {
        double dist = 0;
        for (size_t i = 0; i < dimensions; ++i) {
            double d = (get(i)*pt.get(i)<=1e-20)? 1e10 : log(abs(get(i)/pt.get(i)));
            dist += d * d;
        }
        return dist;
    }
private:
    std::array<coordinate_type, dimensions> coords_;
};
 
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
        double logdistance(const point_type& pt) const {
            return point_.logdistance(pt);
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
 
    node* make_tree(size_t begin, size_t end, size_t index) {
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
 
    void nearest(node* root, const point_type& point, size_t index, bool euc_dist) {
        if (root == nullptr)
            return;
        ++visited_;
        double d = (euc_dist)? root->distance(point): root->logdistance(point);
        if (best_ == nullptr || d < best_dist_) {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->left_ : root->right_, point, index, euc_dist);
        if (dx * dx >= best_dist_)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index, euc_dist);
    }
public:
    kdtree(const kdtree&) = delete;
    kdtree& operator=(const kdtree&) = delete;

	// Default constructor (added by me)
	//kdtree(){};

    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    kdtree(iterator begin, iterator end) : nodes_(begin, end) {
		for ( size_t inode = 0; inode < nodes_.size(); inode++ )
			nodes_[inode].oindex_ = inode; // added by me
        root_ = make_tree(0, nodes_.size(), 0);
    }
  
    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const { return nodes_.empty(); }
 
    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const { return visited_; }
 
    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const { return std::sqrt(best_dist_); }
 
    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& nearest(const point_type& pt, bool euc_dist = true ) {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0, euc_dist);
        return best_->point_;
    }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& nearest(const point_type& pt, size_t & best_oindex, bool euc_dist = true ) {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0, euc_dist);
		best_oindex = best_->oindex_;
        return best_->point_;
    }
};

#endif
