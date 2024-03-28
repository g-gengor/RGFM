#include "RGFM_parser.hpp"

RGFM_parser::RGFM_parser(int grid_num)
{
    this->periodicity_dim = 0;
    this->max_images = Eigen::Vector3i::Zero();
    this->periodic_vecs = Eigen::Matrix3d::Zero();
    this->thread_num = 1;
    this->mpi_partition_num = grid_num;
    this->mode = "INVALID";
    this->gpw_l_max = 8;
        // parse the input file RGFM_config
        ifstream source;
    source.open("RGFM_config", ios_base::in);
    if (!source)
    {
        cerr << "Can't open configuration file";
    }
    else
    {
        for (std::string line; std::getline(source, line);)
        {
            std::istringstream in(line);
            std::string tag;
            in >> tag;
            if (tag == "PERIODIC_DIM")
            {
                in >> this->periodicity_dim;
            }
            else if (tag == "L_MAX")
            {
                in >> this->gpw_l_max;
            }
            else if (tag == "MAX_IMAGES")
            {
                int a, b, c;
                in >> this->max_images(0) >> this->max_images(1) >> this->max_images(2);
                // cout << this->max_images << endl;
            }
            else if (tag == "SPECIES")
            {
                int species_num;
                in >> species_num;
                std::vector<std::string> specs(species_num);
                // cout << species_num << endl;
                for (int i = 0; i < species_num; i++)
                {
                    std::string temp;
                    in >> temp;
                    specs[i] = temp;
                }
                this->gpw_species = specs;
                // cout << this->gpw_species[0] << endl;
                // cout << this->gpw_species[1] << endl;
                // cout << this->gpw_species[2] << endl;
            }
            else if (tag == "PERIODIC_VECS")
            {
                for (int i = 0; i < 3; i++)
                {
                    in >> this->periodic_vecs(i, 0) >> this->periodic_vecs(i, 1) >> this->periodic_vecs(i, 2);
                }
                // cout << this->periodic_vecs << endl;
            }
            else if (tag == "THREAD_NUM")
            {
                in >> this->thread_num;
            }
            else if (tag == "PROBLEM")
            {
                in >> this->mode;
                if (this->mode == "H")
                {
                    cout << "chosen mode: " << this->mode << endl;
                }
                else if (this->mode == "GRID")
                {
                    cout << "chosen mode: " << this->mode << endl;
                }
                else
                {
                    cerr << "Invalid problem identifier";
                }
            }
        }
    }
    //now choose the relevant function and run
    if (this->periodicity_dim == 1 && this->mode == "H")
    {
        //solve for H.txt
        this->disloc_solve_h();
    }
    else if (this->periodicity_dim == 1 && this->mode == "GRID")
    {
        this->disloc_solve_grid();
    }
}

void RGFM_parser::disloc_solve_h()
{
    Greens_PWexpansion G(this->gpw_l_max, this->gpw_species);
    G.set_dislocation_per_vec(this->periodic_vecs(0,0),this->periodic_vecs(0,1),this->periodic_vecs(0,2));//this is the periodicity along dislocation line
    G.set_thread_num(this->thread_num);
    std::cout << "solving for h vector for dislocation" << endl;
    Eigen::VectorXd h = G.solve_h_disloc(this->max_images[0]);
}

void RGFM_parser::disloc_solve_grid()
{
    
    Greens_PWexpansion G(this->gpw_l_max, this->gpw_species);
    G.set_dislocation_per_vec(this->periodic_vecs(0,0),this->periodic_vecs(0,1),this->periodic_vecs(0,2));//this is the periodicity along dislocation line
    G.set_thread_num(this->thread_num);

    std::string grid_prepath("GRIDS/GRID_");
    std::string s = std::to_string(this->mpi_partition_num);
    std::string grid_path = grid_prepath + s;
    std::cout << "Grid solution is starting" << endl;

    auto start = chrono::high_resolution_clock::now();
    Eigen::MatrixX3d u_result = G.solve4grid_disloc_threaded(grid_path,this->max_images[0]);
    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
    cout << "grid solution duration " << duration.count() << " s" << endl;
    std::string U_prepath("US/U_");
    std::string U_path = U_prepath + s;
    std::ofstream file(U_path.c_str());
    if (file.is_open())
    {
        file <<std::setprecision(8)<< u_result << endl;
    }
}