#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

typedef struct Point_
{
    float x, y, z;  // X, Y, Z position
    int clusterID;  // clustered ID
    int original_id; // The id of the range that is related to
}Point;

class DBSCAN {
public:    
    DBSCAN(unsigned int minPts, float eps, const std::vector<Point> &points){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = points;
        m_pointSize = points.size();
    }
    ~DBSCAN(){}

    inline void initialize() {
        for (auto &x:m_points) {
            x.clusterID = UNCLASSIFIED;
        }
    }

    virtual int run();
    inline std::vector<Point> &getPoints() {return m_points;}
    inline int getTotalPointSize() {return m_pointSize;}
    inline int getMinimumClusterSize() {return m_minPoints;}
    inline int getEpsilonSize() {return m_epsilon;}
    inline int getNClusters() {return m_n_clusters;}
    inline std::vector<Point> getCluster(int n_cluster) {
        std::vector<Point> ret;
        if (n_cluster > 0 && n_cluster < m_n_clusters) {
            for (auto &p:m_points) {
                if (p.clusterID == n_cluster) {
                    ret.push_back(p);
                }
            }
        }
        return ret;
    }
protected:
    std::vector<Point> m_points;
    unsigned int m_pointSize;
    unsigned int m_minPoints;
    float m_epsilon;

    virtual std::vector<int> calculateCluster(Point point);
    virtual int expandCluster(Point point, int clusterID);
    inline double calculateDistance(const Point& pointCore, const Point& pointTarget);
    int m_n_clusters;
};

#endif // DBSCAN_H
