/*
 * File: main.cpp
 *
 * @brief A C++ implementation of the the Local Expansion and Conductance
 * Minimizing (LECM) [1]. See the README.md file for compiling and run
 * instructions.
 * @author Guilherme Oliveira Chagas (guilherme.o.chagas[a]gmail.com)
 * @disclaimer Note that I am not a LECM author, so this LECM version may has
 * errors and/or discrepancies with the actual Gao, Zhang and Zhang [1] LECM
 * algorithm.
 * @date This file was created on July 18, 2019, 10:07 AM.
 * @acknowledgment A special thanks to Ph.D. Luis Antonio Nogueira Lorena and
 * Ph.D. Rafael Duarte Coelho dos Santos.
 * @copyright GNU General Public License.
 *
 * References:
 * [1] Y. Gao, H. Zhang and Y. Zhang. Overlapping community detection based on 
 * conductance optimization in large-scale networks, Physica A 522 (2019), 
 * p. 69-79.
 */

#include "../headers/lecm_parameters.hpp"


int main(int argc, char** argv)
{
    Lecm_parameters param;

    if (!param.set_parameters(argc, argv))
    {
        std::cerr << "Please, look at README.md file.\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}