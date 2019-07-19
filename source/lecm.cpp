/*
 * File: lecm.hpp
 * Author: Guilherme O. Chagas
 *
 * @brief Implementation of the Local Expansion and Conductance Minimizing
 * (LECM) [1].
 *
 * (I'm sorry for my bad english xD)
 *
 * Created on May 14, 2019, 04:30 PM
 */

#include "../headers/lecm.hpp"
#include <algorithm>
#include <chrono>
#include <list>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <thread>

// TODO: document this code

/////////////////////////////// Helper functions ///////////////////////////////

namespace
{

/**
* @brief For parallelized methods.
*/
std::mutex my_mutex;

/**
 * @brief
*/
double conductance(const unsigned int g_nb_edges, 
    const unsigned int vol, const unsigned int boundary_edges)
{
    assert(vol <= 2 * g_nb_edges);

    unsigned int compl_vol = 2 * g_nb_edges - vol;
    // prevent division by 0 see LNCS 3418 - Network Analysis
    if (compl_vol == 0 || vol == 0)
    {
        return 1;
    }

    return static_cast<double>(boundary_edges) / std::min(vol, compl_vol);
}

/**
 * @brief
*/
double conductance(const Graph &g, const Cluster &c)
{
    return conductance(g.get_nb_edges(), c.get_volume(), 
        c.get_external_degree());
}


/**
 * @brief
*/
std::pair<unsigned int, unsigned int>
    compute_new_cut_vol(const Cluster &c1, const Cluster &c2, const Graph &g)
{
    unsigned int cut = c2.get_external_degree();
    unsigned int vol = c2.get_volume();
    std::unordered_set<unsigned int> moved_in;

    for (auto v : c1)
    {
        if (!c2.contains(v))
        {
            vol += g.get_vtx_degree(v);
            moved_in.insert(v);
            auto adj_list = g.adj_list_of_vtx(v);
            for (auto it = adj_list.first; it != adj_list.second; ++it)
            {
                if (c2.contains(*it))
                {
                    --cut; // then, decrease the number of boundary edges
                }
                else if (moved_in.find(*it) != moved_in.end()) // was moved in
                {
                    --cut; // then the number of boundary edges is decreased
                }
                else // the adj vtx is not contained in cluster
                {
                    ++cut; // then, increase the number of boundary edges
                }
            }
        }
    }

    return std::make_pair(cut, vol);
}


/**
 * @brief Check whether the conductance of the cluster decreased or not after 
 * the removal of the vertex. The cluster cut and volume are passed by reference
 * and are updated inside this fuction if a conductance decreased is detected.
 * This method is used in the node movement step (see [1,2]).
 * @param unsigned int &: cluster's currently cut (external degree).
 * @param unsigned int &: cluster's currently volume.
 * @param const unsigned int: vertex id to be removed.
 * @param const std::unordered_set<unsigned int> &: set of vertices marked as 
 * "move out", that is, vertices that are still in the cluster but in practice 
 * were removed from it.
 * @param const Cluster &: cluster c.
 * @param const Graph&: graph g.
 * @return bool: true with the cluster's conductance was decreased after v 
 * removal. False, otherwise.
*/
bool conductance_r_decreased(unsigned int &cut, unsigned int &vol, 
    const unsigned int v, const std::unordered_set<unsigned int> &move_out,
    const Cluster &c, const Graph& g)
{
    unsigned int new_cut = cut; // new number of boundary edges
    auto adj_list = g.adj_list_of_vtx(v);
    for (auto it = adj_list.first; it != adj_list.second; ++it)
    {
        if (c.contains(*it))
        {
            // if adj vtx is contained in the cluster and was "moved out"
            if (move_out.find(*it) != move_out.end())
            {
                --new_cut; // then reduce the number of boundary edges
            }
            else // otherwise the adj vtx is contained and not moved out yet
            {
                ++new_cut;
            }
        }
        else // if adj vtx is not contained
        {
            --new_cut; // then edge is not a boundary anymore
        }
    }

    double orig_cond = conductance(g.get_nb_edges(), vol, cut);
    double new_cond = conductance(g.get_nb_edges(), 
        vol - g.get_vtx_degree(v), new_cut);

    if ((orig_cond - new_cond) > 0) // if conductance was decreased
    {
        cut = new_cut; // update cluster cut after removing v
        vol -= g.get_vtx_degree(v); // update cluster volume after removing v
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * @brief Check whether the conductance of the cluster decreased or not after 
 * the insertion of the vertex. The cluster cut and volume are passed by 
 * reference and are updated inside this fuction if a conductance decreased is 
 * detected. This method is used in the node movement step (see [1,2]).
 * @param unsigned int &: cluster's currently cut (external degree).
 * @param unsigned int &: cluster's currently volume.
 * @param const unsigned int: vertex id to be inserted.
 * @param const std::unordered_set<unsigned int> &: set of vertices marked as 
 * "move in", that is, vertices that are not contained in the cluster but in 
 * practice they are.
 * @param const std::unordered_set<unsigned int> &: set of vertices marked as 
 * "move out", that is, vertices that are still in the cluster but in practice 
 * were removed from it.
 * @param const Cluster &: cluster c.
 * @param const Graph&: graph g.
 * @return bool: true with the cluster's conductance was decreased after v 
 * insertion. False, otherwise.
*/
bool conductance_i_decreased(unsigned int &cut, unsigned int &vol,
    const unsigned int v, const std::unordered_set<unsigned int> &move_in, 
    const std::unordered_set<unsigned int> &move_out, const Cluster &c, 
    const Graph& g)
{
    unsigned int new_cut = cut; // new number of boundary edges
    auto adj_list = g.adj_list_of_vtx(v);
    for (auto it = adj_list.first; it != adj_list.second; ++it)
    {
        if (c.contains(*it))
        {
            // if adj vtx is contained in the cluster and was "moved out"
            if (move_out.find(*it) != move_out.end()) 
            {
                ++new_cut; // then increase the number of boundary edges
            }
            else // otherwise, the vtx is indeed contained
            {
                --new_cut; // then, decrease the number of boundary edges
            }
        }
        else if (move_in.find(*it) != move_in.end()) // adj vtx was moved in
        {
            --new_cut; // then the number of boundary edges is decreased
        }
        else // the adj vtx is not contained in cluster
        {
            ++new_cut; // then, increase the number of boundary edges
        }
    }

    double orig_cond = conductance(g.get_nb_edges(), vol, cut);
    double new_cond = conductance(g.get_nb_edges(), 
        vol + g.get_vtx_degree(v), new_cut);

    if ((orig_cond - new_cond) > 0) // if conductance was decreased
    {
        cut = new_cut; // update cluster cut after removing v
        vol += g.get_vtx_degree(v); // update cluster volume after removing v
        return true;
    }
    else
    {
        return false;
    }
}

/** 
 * @brief Equation 10 of [1].
 * @param
 * @param
 * @return
*/
unsigned int equation_10(
    const std::unordered_map<unsigned int, unsigned int> clsts_nb_edges,
    const unsigned int v_dg)
{
    unsigned c_id = 0;
    double max_ratio = std::numeric_limits<double>::lowest();
    for (auto itc = clsts_nb_edges.begin(); itc != clsts_nb_edges.end(); ++itc)
    {
        double ratio = static_cast<double>((*itc).second) / v_dg;
        if (ratio > max_ratio)
        {
            max_ratio = ratio;
            c_id = (*itc).first;
        }
    }
    return c_id;
}

/** 
 * @brief: map of the adjacent clsts and the number of edges between them key is
 * the cluster index and the mapped value is the number of edges between the 
 * vertex v and the cluster.
 * @param:.
 * @param:.
 * @param:.
 * @return:.
*/
std::unordered_map<unsigned int, unsigned int> edges_between_v_and_clsts(
    const Graph &g,
    const unsigned int v,
    const Clustering &clusters)
{
    std::unordered_map<unsigned int, unsigned int> edges_v_c;
    auto adj_list = g.adj_list_of_vtx(v);
    for (auto ita = adj_list.first; ita != adj_list.second; ++ita)
    {
        auto v_belonging_c = clusters.get_v_belonging(*ita);
        for (auto itc = v_belonging_c.first; itc != v_belonging_c.second; ++itc)
        {
            if (edges_v_c.find(*itc) == edges_v_c.end())
            {
                edges_v_c.insert({*itc, 1});
            }
            else
            {
                ++edges_v_c.at(*itc);
            }
        }
    }
    return edges_v_c;
}

/** 
 * @brief Inequality 9 of [1].
 * @param
 * @param
 * @param
 * @return
*/
bool inequality_9(const double chi, const unsigned int edges_v_ci, 
    const unsigned int v_dg)
{
    return (static_cast<double>(edges_v_ci) / v_dg) > chi;
}

/** 
 * @brief Inequality 11 of [1].
 * @param
 * @param
 * @param
 * @return
*/
bool inequality_11(const Cluster &c, const unsigned int v_dg, 
    const unsigned int edges_v_ci)
{
    double lhs = 1 - 2 * static_cast<double>(edges_v_ci) / v_dg;

    double rhs_dividend = c.get_volume() - 2 * c.get_internal_degree();
    double rhs = rhs_dividend / c.get_volume();

    return lhs < rhs;
}

/**
* @brief Sort vertices indices in decreasing probability-per-degree (PPD) order.
* This fuction returns a vector of vertices indices sorted in decreasing 
* probability-per-degree. See [1,2] for details.
* @param: std::unordered_map<unsigned int, double> &.
* @param: Graph &:.
*/
std::vector<unsigned int> sort_vertices_in_decreasing_ppd(
    const std::unordered_map<unsigned int, double> &p, const Graph &g)
{
    std::vector<unsigned int> decreasing_ppd(p.size()); // vertices indices
    unsigned int i = 0;
    for (auto &e : p) // populate vector with vertices indices
    {
        decreasing_ppd[i++] = e.first;
    }
    std::sort(decreasing_ppd.begin(), decreasing_ppd.end(), // sort
        [&](const unsigned int v1, const unsigned int v2)
        {
            return (p.at(v1) / static_cast<double>(g.get_vtx_degree(v1))) 
                > (p.at(v2) / static_cast<double>(g.get_vtx_degree(v2)));
        });
    return decreasing_ppd;
}

/**
* @brief Sort vertices indices in decreasing probability-per-degree (PPD) order.
* This fuction returns a vector of vertices indices sorted in decreasing 
* probability-per-degree. See [1,2] for details.
* @param: std::unordered_map<unsigned int, double> &.
* @param: Graph &:.
*/
Cluster sweep_step(const Graph &g, const std::vector<unsigned int> &ppd_order)
{
    Cluster clst(g); // cluster constructed from original graph
    // std::set<unsigned int> sweep_set;
    std::unordered_map<unsigned int, bool> swept;
    unsigned int vol = 0;
    unsigned int cut = 0; // cluster boundary edges
    double min_cond = std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < ppd_order.size(); ++i)
    {
        // compute the number of boundary edges
        swept[ppd_order[i]] = true;
        auto adj_list = g.adj_list_of_vtx(ppd_order[i]);
        for (auto it_adj = adj_list.first; it_adj != adj_list.second; ++it_adj)
        {
            if (swept[*it_adj])
            {
                --cut;
            }
            else
            {
                ++cut;
            }
        }
        vol += g.get_vtx_degree(ppd_order[i]);
        // compute conductance
        double cond = conductance(g.get_nb_edges(), vol, cut);
        if (cond <= min_cond)
        {
            min_cond = cond;
            // save set of vertices (clst) of minimum conductance
            for (unsigned int j = clst.size(); j <= i; ++j)
            {
                // insert vertices from the original graph
                clst.insert(ppd_order[j]);
            }
        }
    }

    return clst;
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////


Lecm::Lecm(const Lecm_parameters &p, const Graph &g) :
    m_p(p),
    m_graph(g),
    m_executed(false),
    m_clusters(g.get_nb_vertices())
{}


Lecm::~Lecm() {}


void Lecm::execute()
{
    std::cout << "Executing LECM...\n";
    auto start = std::chrono::high_resolution_clock::now();
    auto seeds = seeding();

    seed_expansion(seeds);

    community_refinement();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "It took me " << 
        std::chrono::duration_cast<std::chrono::duration<double> >
        (end - start).count() << "s.\n";
    m_executed = true;
}


void Lecm::write_clustering() const
{
    std::cout << "Writing clustering in file...\n";
    if (!m_executed)
    {
        std::cerr << "LECM algorithm hasn't been executed! I can't write any "
            "clustering! Aborting...\n";
        exit(EXIT_FAILURE);
    }

    m_clusters.write_clusters_in_file("./clustering.dat");
}


/////////////////////////////// private methods ////////////////////////////////


void Lecm::community_refinement()
{
    std::cout << "\tCommunity refinement step...\n";
    node_movement();

    combine_communities();

    find_communities_for_outliers();
}


double Lecm::cs(const unsigned int c_i, const unsigned int c_j) const
{
    unsigned int cut = 0;
    unsigned int vol = 0;

    const Cluster &c1 = m_clusters[c_i];
    const Cluster &c2 = m_clusters[c_j];

    if (c1.size() < c2.size())
    {
        std::tie(cut, vol) = compute_new_cut_vol(c1, c2, m_graph);
    }
    else
    {
        std::tie(cut, vol) = compute_new_cut_vol(c1, c2, m_graph);
    }

    double sum = conductance(m_graph, c1) + conductance(m_graph, c2);
    double newc = conductance(m_graph.get_nb_edges(), vol, cut);

    unsigned int inter_size = m_clusters.get_inter_size(c_i, c_j);
    unsigned int union_size = c1.size() + c2.size() - inter_size;
    double jaccard = static_cast<double>(inter_size) / union_size;

    return (m_p.theta() * jaccard) + (1 - m_p.theta()) * sum / (2 * newc + sum);
}


void Lecm::node_movement()
{
    #pragma omp parallel for
    for (unsigned int i = 0; i < m_clusters.size(); ++i)
    {
        unsigned int vol = m_clusters[i].get_volume();
        unsigned int cut = m_clusters[i].get_external_degree();

        // find vertices to move out (see [2])
        std::unordered_set<unsigned int> move_out;
        for (auto v : m_clusters[i])
        {
            // if true, cut and vol are updated inside "conductance_r_decreased"
            if (conductance_r_decreased(cut, vol, v, move_out, m_clusters[i],
                m_graph))
            {
                move_out.insert(v);
            }
        }

        // find vertices to move in
        std::unordered_set<unsigned int> move_in;
        auto ext_adj_v = m_clusters[i].get_external_adj_vertices();
        for (auto it = ext_adj_v.first; it != ext_adj_v.second; ++it)
        {
            if (conductance_i_decreased(cut, vol, *it, move_in, move_out, 
                m_clusters[i], m_graph))
            {
                move_in.insert(*it);
            }
        }

        // refine cluster C moving vertices out and in
        for (auto v : move_out)
        {
            #pragma omp critical
            {
                m_clusters.remove_v_from_clst(v, i);
            }
        }

        for (auto v : move_in)
        {
            #pragma omp critical
            {
                m_clusters.insert_v_in_clst(v, i);
            }
        }
    }
}


void Lecm::combine_communities()
{
    std::unordered_set<unsigned int> deleted_clusters; // workaround
    // iterate over the Com_neighbors list until the current last cluster
    for (auto it = m_clusters.begin(); it != m_clusters.end(); ++it)
    {
        unsigned int c_i = (*it).first;
        if (deleted_clusters.find(c_i) != deleted_clusters.end())
        {
            continue; // "empty/inserted" cluster, get next iterator
        }

        auto oc = m_clusters.get_ovlp_clsts(c_i); // ovlp clusters
        // for each cluster that overlaps
        for (auto itoc = oc.first; itoc != oc.second; )
        {
            unsigned int c_j = (*itoc).first;
            // checks if an improvement is found by merging clusters
            if (deleted_clusters.find(c_j) == deleted_clusters.end() &&
                cs(c_i, c_j) >= m_p.beta())
            {
                // an improvement is detected then delete c_j is delete and
                // iterates over c_j vertices removing they and inserting in c_i
                deleted_clusters.insert(c_j);
                for (auto itv = m_clusters[c_j].begin(); 
                    itv != m_clusters[c_j].end(); )
                {
                    m_clusters.insert_v_in_clst(*itv, c_i);
                    itv = m_clusters.remove_v_from_clst(itv, c_j); // next it
                }
                // get new overlapping clusters list
                oc = m_clusters.get_ovlp_clsts(c_i);
                itoc = oc.first; // starts from begin again
            }
            else // if an improvement is not found, compare next ovlp cluster
            {
                ++itoc;
            }
        }
    }

    /* remove the "empty/inserted" clusters. Merged clusters are not removed 
    when they are merged because we are iterating over the m_clusters contained
    and a removal may invalidad the its iterators. So, we keep an "empty shell"
    and remove it they now when is safe.*/
    for (auto c : deleted_clusters)
    {
        m_clusters.remove_clst(c);
    }
}


void Lecm::find_communities_for_outliers()
{
    std::vector<unsigned int> remaning_v; // vertices that were not inserted
    for (unsigned int v = 0; v < m_graph.get_nb_vertices(); ++v)
    {
        if (m_clusters.get_v_belonging_size(v) == 0)
        {
            // get the edges between vertex v and all clusters
            auto edges_v_c = edges_between_v_and_clsts(m_graph, v, m_clusters);

            for (auto itc = edges_v_c.begin(); itc != edges_v_c.end(); ++itc)
            {
                const Cluster &c = m_clusters[(*itc).first];
                if (inequality_11(c, m_graph.get_vtx_degree(v), (*itc).second))
                {
                    m_clusters.insert_v_in_clst(v, (*itc).first);
                }
            }

            // if still not belongs to none
            if (m_clusters.get_v_belonging_size(v) == 0)
            {
                remaning_v.push_back(v);
            }
        }
    }

    // try to find a cluster for the vertices that still not belong to none
    for (auto it = remaning_v.begin(); it != remaning_v.end(); ++it)
    {
        auto edges_v_c = edges_between_v_and_clsts(m_graph, *it, m_clusters);

        if (edges_v_c.size() > 0)
        {
            unsigned int v_dg = m_graph.get_vtx_degree(*it);
            unsigned int c_id = equation_10(edges_v_c, v_dg);
            if (inequality_9(m_p.chi(), edges_v_c.at(c_id), v_dg))
            {
                m_clusters.insert_v_in_clst(*it, c_id);
            }
        }
    }
}


void Lecm::seed_expansion(const std::vector<std::vector<unsigned int>> &seeds)
{
    std::cout << "\tSeed expansion step...\n";
    const unsigned int nb_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.reserve(nb_threads);
    m_clusters.reserve(seeds.size());

    for (unsigned int i = 0; i < nb_threads; ++i)
    {
        threads.push_back(std::thread(&Lecm::seed_expansion_thread_task, this,
            std::ref(seeds), i));
    }

    for (auto &t : threads)
    {
        t.join();
    }
}


void Lecm::seed_expansion_thread_task(
    const std::vector<std::vector<unsigned int>> &seeds,
    const unsigned int thread_id)
{
    const unsigned int nb_threads = std::thread::hardware_concurrency();
    const unsigned int thread_offset = seeds.size() / nb_threads;
    const unsigned int rem = seeds.size() % nb_threads;
    unsigned int starting_seed = thread_id * thread_offset; 
    unsigned int final_seed = starting_seed + thread_offset;

    if (thread_id >= nb_threads - rem) // remaining seeds
    {
        starting_seed += thread_id - (nb_threads - rem);
        final_seed += thread_id - (nb_threads - rem - 1);
    }

    for (unsigned int i = starting_seed; i < final_seed; ++i)
    {
        auto p = ppr(seeds[i]);
        auto ppd_order = sort_vertices_in_decreasing_ppd(p, m_graph);
        Cluster clst(sweep_step(m_graph, ppd_order));

        if (clst.size() > 2) // keep only clusters with 3 or more vertices [1]
        {
            std::lock_guard<std::mutex> locker(my_mutex);
            m_clusters.insert_without_repetition(clst);
        }
    }
}


std::vector<std::vector<unsigned int>> Lecm::seeding() const
{
    std::cout << "\tSeeding step...\n";
    std::vector<std::vector<unsigned int>> seeds;
    std::vector<bool> taken(m_graph.get_nb_vertices(), false);
    std::vector<unsigned int> order(m_graph.get_nb_vertices());

    // populate order array with vertices id (from 0 to n-1)
    std::iota(order.begin(), order.end(), 0);
    // sort vertices decreasingly by their degree
    std::sort(order.begin(), order.end(), 
        [&](const unsigned int v, const unsigned int u)
        {
            return m_graph.get_vtx_degree(v) > m_graph.get_vtx_degree(u);
        });

    for (auto v : order)
    {
        if (!taken[v])
        {
            taken[v] = true;
            std::vector<unsigned int> seed;
            seed.reserve(m_graph.get_vtx_degree(v) + 1);
            seed.push_back(v);
            // mark neighboors as taken and insert them in seed vector
            auto adj_list = m_graph.adj_list_of_vtx(v);
            for (auto it = adj_list.first; it != adj_list.second; ++it)
            {
                if (!taken[*it])
                {
                    taken[*it] = true;
                    seed.push_back(*it);
                }
            }
            seeds.push_back(seed);
        }
    }

    return seeds;
}


std::unordered_map<unsigned int, double> 
    Lecm::ppr(const std::vector<unsigned int> &seed) const
{
    std::unordered_map<unsigned int, double> p;
    std::unordered_map<unsigned int, double> r;

    std::queue<unsigned int> q; // queue of vertices as presented by [3]
    // initialize p and r values and stores vertices in the queue
    for (auto v : seed)
    {
        p[v] = 0;
        r[v] = static_cast<double>(1) / seed.size();
        // if r[v] > degree(v) * epsilon, then push to queue
        if (r[v] > m_graph.get_vtx_degree(v) * m_p.epsilon())
        {
            q.push(v);
        }
    }

    while (!q.empty()) // compute and update x and r values
    {
        unsigned int v = q.front();
        p[v] += (1 - m_p.alpha()) * r[v];
        auto adj_list = m_graph.adj_list_of_vtx(v);
        for (auto it = adj_list.first; it != adj_list.second; ++it)
        {
            r[*it] += m_p.alpha() * r[v] / (2 * m_graph.get_vtx_degree(v));
            if (r[*it] > m_graph.get_vtx_degree(*it) * m_p.epsilon())
            {
                q.push(*it);
            }
        }
        r[v] = m_p.alpha() * r[v] / 2;
        // remove from queue front vertices with r[v] <= degree(v) * epsilon
        while (r[q.front()] <= m_graph.get_vtx_degree(q.front()) *m_p.epsilon())
        {
            q.pop();
            if (q.empty())
            {
                break;
            }
        }
    }

    // sweep step: pick and return the set (cluster) with minimum conductance
    return p;
}