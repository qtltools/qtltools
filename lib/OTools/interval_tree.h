/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#ifndef _INTERVAL_TREE_H
#define _INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;


template <class T, typename K = int>
class Interval {
public:
    K start;
    K stop;
    T value;
    Interval(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};

template <class T, typename K>
int intervalStart(const Interval<T,K>& i) {
    return i.start;
}

template <class T, typename K>
int intervalStop(const Interval<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
ostream& operator<<(ostream& out, Interval<T,K>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = int>
class IntervalStartSorter {
public:
    bool operator() (const Interval<T,K>& a, const Interval<T,K>& b) {
        return a.start < b.start;
    }
};

template <class T, typename K = int>
class IntervalTree {

public:
    typedef Interval<T,K> interval;
    typedef vector<interval> intervalVector;
    typedef IntervalTree<T,K> intervalTree;

    intervalVector intervals;
    intervalTree* left;
    intervalTree* right;
    int center;

    IntervalTree(void)
        : left(NULL)
        , right(NULL)
        , center(0)
    { }

    IntervalTree(const intervalTree& other) {
        center = other.center;
        intervals = other.intervals;
        if (other.left) {
            left = (intervalTree*) malloc(sizeof(intervalTree));
            *left = *other.left;
        } else {
            left = NULL;
        }
        if (other.right) {
            right = new intervalTree();
            *right = *other.right;
        } else {
            right = NULL;
        }
    }

    IntervalTree(
            intervalVector& ivals,
            unsigned int depth = 16,
            unsigned int minbucket = 64,
            int leftextent = 0,
            int rightextent = 0,
            unsigned int maxbucket = 512
            )
        : left(NULL)
        , right(NULL)
    {

        --depth;
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            intervals = ivals;
        } else {
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
                IntervalStartSorter<T,K> intervalStartSorter;
                sort(ivals.begin(), ivals.end(), intervalStartSorter);
            }

            int leftp = 0;
            int rightp = 0;
            int centerp = 0;

            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                vector<K> stops;
                stops.resize(ivals.size());
                transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
                rightp = *max_element(stops.begin(), stops.end());
            }

            //centerp = ( leftp + rightp ) / 2;
            centerp = ivals.at(ivals.size() / 2).start;
            center = centerp;

            intervalVector lefts;
            intervalVector rights;

            for (typename intervalVector::iterator i = ivals.begin(); i != ivals.end(); ++i) {
                interval& interval = *i;
                if (interval.stop < center) {
                    lefts.push_back(interval);
                } else if (interval.start > center) {
                    rights.push_back(interval);
                } else {
                    intervals.push_back(interval);
                }
            }

            if (!lefts.empty()) {
                left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
            }
            if (!rights.empty()) {
                right = new intervalTree(rights, depth, minbucket, centerp, rightp);
            }
        }
    }

    IntervalTree& operator=(const intervalTree& other) {
        center = other.center;
        intervals = other.intervals;
        if (other.left) {
            left = new intervalTree();
            *left = *other.left;
        } else {
            left = NULL;
        }
        if (other.right) {
            right = new intervalTree();
            *right = *other.right;
        } else {
            right = NULL;
        }
        return *this;
    }

    void findOverlapping(K start, K stop, intervalVector& overlapping) {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
                interval& interval = *i;
                if (interval.stop >= start && interval.start <= stop) {
                    overlapping.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findOverlapping(start, stop, overlapping);
        }

        if (right && stop >= center) {
            right->findOverlapping(start, stop, overlapping);
        }

    }

    void findOverlapping(K pos, intervalVector& overlapping) {
        if (!intervals.empty() && ! (pos < intervals.front().start)) {
            for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
                interval& interval = *i;
                if (interval.stop >= pos && interval.start <= pos) {
                    overlapping.push_back(interval);
                }
            }
        }

        if (left && pos <= center) {
            left->findOverlapping(pos, overlapping);
        }

        if (right && pos >= center) {
            right->findOverlapping(pos, overlapping);
        }

    }

    void findContained(K start, K stop, intervalVector& contained) {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
                interval& interval = *i;
                if (interval.start >= start && interval.stop <= stop) {
                    contained.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findContained(start, stop, contained);
        }

        if (right && stop >= center) {
            right->findContained(start, stop, contained);
        }

    }

    bool checkOverlapping(K pos) {
    	bool outcome = false;

    	if (!intervals.empty() && ! (pos < intervals.front().start)) {
    		for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
    			if (i->stop >= pos && i->start <= pos)
    				outcome = true;
            }
    	}

        if (left && pos <= center) {
        	outcome = outcome || left->checkOverlapping(pos);
        }

        if (right && pos >= center) {
        	outcome = outcome || right->checkOverlapping(pos);
        }
        return outcome;
    }

    ~IntervalTree(void) {
        // traverse the left and right
        // delete them all the way down
        if (left) {
            delete left;
        }
        if (right) {
            delete right;
        }
    }

};

#endif
