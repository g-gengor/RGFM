#include "Greens_PWexpansion.hpp"

#define MAXBUFSIZE  ((int) 1e6)

Greens_PWexpansion::Greens_PWexpansion(int l_max, std::vector<std::string> species)
{
    this->thread_num_set = false;
    
    this->gpw_l_max = l_max;
    this->gpw_species = species;
    this->gpw_total_pw_loop_size = ((gpw_l_max / 2) + 1) * ((gpw_l_max / 2) + 1);
    //read the anisotropy terms
    cout << "reading A" << endl;
    string A_path_append_real("_r.txt");
    string A_path_append_imag("_i.txt");
    this->A_matrices.resize(this->gpw_total_pw_loop_size * 3,3);
    int m_order_count = 0;
    for (int l = 0; l < gpw_l_max + 1; l += 2)
    {
        string l_num(std::to_string(l));
        for (int m = 0; m <= l; m++)
        {
            string A_path_real("A_matrices/A_");
            string A_path_imag("A_matrices/A_");
            A_path_real.append(l_num);
            A_path_real.append("_");
            A_path_imag.append(l_num);
            A_path_imag.append("_");
            string m_num(std::to_string(m));
            A_path_real.append(m_num);
            A_path_real.append(A_path_append_real);
            A_path_imag.append(m_num);
            A_path_imag.append(A_path_append_imag);

            MatrixXd realpart = readMatrix(A_path_real.c_str());
            MatrixXd imagpart = readMatrix(A_path_imag.c_str());

            this->A_matrices.block<3,3>(m_order_count*3,0) = realpart + imagpart * 1i;
            m_order_count++;
        }
        
    }
    cout << "A is read" << endl;

    //read the force profile data
    cout << "reading F" << endl;
    int num_of_species = this->gpw_species.size();
    for (int i = 0; i < num_of_species; i++)
    {
        std::string sp_name = this->gpw_species[i];
        std::string force_path("F/F_");
        force_path.append(sp_name);
        force_path.append(".txt");
        MatrixXd force_profile = readMatrix(force_path.c_str());
        this->F_l.push_back(force_profile);
    }
    cout << "F read successfully" << endl;
    
    //read the atomic displacement data
    cout << "reading DISP" << endl;
    MatrixXd disps_raw = readMatrix("DISP");
    this->disps.resize(disps_raw.rows(),7);
    MatrixXd temp = disps_raw.block(0,1,disps_raw.rows(),7);
    MatrixXd temp2 = temp;
    temp2.block(0,1,temp.rows(),3) = temp.block(0,1,temp.rows(),3) - temp.block(0,4,temp.rows(),3);
    this->disps = temp2;
    this->total_atom_num = this->disps.rows();
    cout << "total atom number: " << this->total_atom_num << endl;
    cout << "DSIP is read" << endl;
    cout << "DISP: " << this->disps.rows() << "x" << this->disps.cols() << endl;

    //std::ofstream file("A_r.txt");
    //if (file.is_open())
    //{
    //    file << this->A_matrices.real() << endl;
    //}
    //std::ofstream file2("A_i.txt");
    //if (file2.is_open())
    //{
    //    file2 << this->A_matrices.imag() << endl;
    //}
}

Greens_PWexpansion::Greens_PWexpansion(int l_max, std::vector<std::string> species, bool system_solved)
{
    this->thread_num_set = false;
    
    this->gpw_l_max = l_max;
    this->gpw_species = species;
    this->gpw_total_pw_loop_size = ((gpw_l_max / 2) + 1) * ((gpw_l_max / 2) + 1);
    //read the anisotropy terms
    cout << "reading A" << endl;
    string A_path_append_real("_r.txt");
    string A_path_append_imag("_i.txt");
    this->A_matrices.resize(this->gpw_total_pw_loop_size * 3,3);
    int m_order_count = 0;
    for (int l = 0; l < gpw_l_max + 1; l += 2)
    {
        string l_num(std::to_string(l));
        for (int m = 0; m <= l; m++)
        {
            string A_path_real("A_matrices/A_");
            string A_path_imag("A_matrices/A_");
            A_path_real.append(l_num);
            A_path_real.append("_");
            A_path_imag.append(l_num);
            A_path_imag.append("_");
            string m_num(std::to_string(m));
            A_path_real.append(m_num);
            A_path_real.append(A_path_append_real);
            A_path_imag.append(m_num);
            A_path_imag.append(A_path_append_imag);

            MatrixXd realpart = readMatrix(A_path_real.c_str());
            MatrixXd imagpart = readMatrix(A_path_imag.c_str());

            this->A_matrices.block<3,3>(m_order_count*3,0) = realpart + imagpart * 1i;
            m_order_count++;
        }
        
    }
    cout << "A is read" << endl;

    //read the force profile data
    cout << "reading F" << endl;
    int num_of_species = this->gpw_species.size();
    for (int i = 0; i < num_of_species; i++)
    {
        std::string sp_name = this->gpw_species[i];
        std::string force_path("F/F_");
        force_path.append(sp_name);
        force_path.append(".txt");
        MatrixXd force_profile = readMatrix(force_path.c_str());
        this->F_l.push_back(force_profile);
    }
    cout << "F are read" << endl;
    
    //read the atomic displacement data
    cout << "reading DISP" << endl;
    MatrixXd disps_raw = readMatrix("DISP");
    this->disps.resize(disps_raw.rows(),7);
    MatrixXd temp = disps_raw.block(0,1,disps_raw.rows(),7);
    MatrixXd temp2 = temp;
    temp2.block(0,1,temp.rows(),3) = temp.block(0,1,temp.rows(),3) - temp.block(0,4,temp.rows(),3);
    this->disps = temp2;
    this->total_atom_num = this->disps.rows();
    cout << "DSIP is read" << endl;
    cout << "DISP: " << this->disps.rows() << "x" << this->disps.cols() << endl;

    //read the fitting coefficients
    cout << "reading H" << endl;
    this->h = readMatrix("H.txt");
    cout << "H is read" << endl;
    cout << "H: " << this->h.rows() << "x" << this->h.cols() << endl;
}

Greens_PWexpansion::~Greens_PWexpansion()
{
}

MatrixXd Greens_PWexpansion::readMatrix(const char *filename)
    {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
        {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
    }

Matrix3d Greens_PWexpansion::getGreens(int atom_index, double px, double py, double pz){
    //intiliaze a complex amtrix for computations
    Matrix3cd Greens_temp = Matrix3cd::Zero();
    //get the position atom that we are calculation greens function for
    Vector3d atom_pos;
    atom_pos(0) = this->disps(atom_index,1);
    atom_pos(1) = this->disps(atom_index,2);
    atom_pos(2) = this->disps(atom_index,3);
    //cout << atom_pos << endl;
    //atom_pos.transposeInPlace();
    //calcualte the vector between enquired point and atom source
    Vector3d position;
    position(0) = px;
    position(1) = py;
    position(2) = pz;
    Vector3d x = position - atom_pos;
    //cout << atom_pos << endl;
    int species = this->disps(atom_index,0) - 1;//get the species enumerator
    //cout << "species: " << species << endl;
    double r = x.norm();
    double theta, phi;
    //cout << "x: " << x << endl;
    //cout << "r: " << r << endl;
    if (r == 0.0)
    {
        theta = 0.0;
        phi = 0.0;
    }
    else
    {
        theta = std::acos(x(2) / r);
        phi = atan2(x(1), x(0));
        if (phi < 0)
        {
            phi = phi + 2 * M_PI;
        }
        if (theta < 0)
        {
            theta = theta + M_PI;
        }
        
        
    }
    //cout << "t: " << theta << endl;
    //cout << "p: " << phi << endl;
    //now start the loop for calculation
    //following is for the linear interppolation of F_ls
    int m_counter = 0;
    //check the lower and upper bounds for F_l
    double lb_F = F_l[species](0,0);
    double ub_F = F_l[species].col(0).tail<1>().value();
    bool F_inbound = 0;
    if (r >= lb_F && r <= ub_F)
    {
        F_inbound = 1;
    }
    //cout << "G: " << endl;
    //cout << Greens_temp << endl;
    for (int l = 0; l <= this->gpw_l_max; l += 2)
    {
        //l goes 0,2,4,6,8...
        //calculate the terms that are only dependent on l
        std::complex<double> c_temp = 1i;
        std::complex<double> c = pow(c_temp,l * 1.0); //1.0 is there to make sure l is converted to double
        //interpoalte F_l. currently interpoaltion order is 1
        double F = 0.0;
        if (F_inbound)
        {
            std::vector<double> F_l_x;
            std::vector<double> F_l_y;
            //resize the holder vectors
            F_l_x.resize(F_l[species].col(0).size());
            F_l_y.resize(F_l[species].col(0).size());
            //assign the values to the holder vectors
            VectorXd::Map(&F_l_y[0], F_l[species].col((l / 2) + 1).size()) = F_l[species].col((l / 2) + 1);
            VectorXd::Map(&F_l_x[0], F_l[species].col(0).size()) = F_l[species].col(0);//first column is the x values
            //boost interpolation library
            boost::math::interpolators::barycentric_rational<double> interpolant(std::move(F_l_x), std::move(F_l_y), 1 );
            F = interpolant(r);
        }
        else if (r < lb_F)
        {
            F = F_l[species].col((l / 2) + 1)(0);
        }
        else
        {
            F = 0;//F_l[species].col((l / 2) + 1).tail<1>().value();
        }        
        //cout << F << endl;
        //cout << c << endl;
        for (int m = 0; m <= l; m++)
        {
            double m_factor = 1.0;
            if (m > 0)
            {
                // when m is not zero -m and m yields the same result for Greens tensor
                m_factor = 2.0;
            }
            
            std::complex<double> Ylm = boost::math::spherical_harmonic(l,m,theta,phi);
            //cout << Ylm << endl;
            //get the corresponding Anisotropy matrix
            Matrix3cd Alm = this->A_matrices.block(m_counter * 3, 0, 3, 3);
            Matrix3cd placeholder = m_factor * c * Ylm * F * Alm / (2 * M_PI * M_PI) + Greens_temp;
            Greens_temp = placeholder;
            //cout << "c: " << c << endl;
            //cout << "Ylm: " << Ylm << endl;
            //cout << "F: " << F << endl;
            //cout << "Alm" << endl;
            //cout << Alm << endl;
            //cout << "m factor: " << m_factor << endl;
            //cout << "G: " << endl;
            //cout << Greens_temp << endl;
            m_counter++;
        }
    }
    return Greens_temp.real();
}

VectorXd Greens_PWexpansion::solve_h(){
    //this function finds an h vector. However, when it writes it, it is just bogus. I could not solve for the life of me
    cout << "beginning of solve_h" << endl;
    //get the displacements we are interested in
    VectorXd u(this->total_atom_num*3);
    for (int i = 0; i < this->total_atom_num; i++)
    {
        u(i*3) = this->disps(i,4);
        u(i*3+1) = this->disps(i,5);
        u(i*3+2) = this->disps(i,6);
    }
    
    /*std::ofstream file3("solve4.txt");
    if (file3.is_open())
    {
        file3 << u << endl;
    }*/
    //cout << "original:" << endl << this->disps.block(0,4,this->total_atom_num,3) << endl;
    //cout << "reshaped:" << endl << this->disps.block(0,4,this->total_atom_num,3).reshaped<RowMajor>() << endl;
    //cout << u << endl;

    //build the system matrix
    MatrixXd system_matrix(this->total_atom_num*3, this->total_atom_num*3);
    cout << "system matrix loop begins" << endl;
    auto start = chrono::high_resolution_clock::now();
    for (int t_ind = 0; t_ind < this->total_atom_num; t_ind++)
    {
        Vector3d t_pos;
        t_pos(0) = this->disps(t_ind,1);
        t_pos(1) = this->disps(t_ind,2);
        t_pos(2) = this->disps(t_ind,3);
        for (int s_ind = 0; s_ind < this->total_atom_num; s_ind++)
        {
            Matrix3d G_partial = this->getGreens(s_ind,t_pos(0), t_pos(1), t_pos(2));
            system_matrix.block(t_ind*3,s_ind*3,3,3) = G_partial;
        }
        
    }
    
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "system matrix is created in " << duration.count() << " ms" << endl;
    cout << "system matrix shape rows x cols: " << system_matrix.rows() << " x " << system_matrix.cols() << endl;
    cout << "u shape: " << u.rows() << " x " << u.cols() << endl;
    //solve the linear system for the fitting coefficients
    
    Eigen::FullPivLU<MatrixXd> fpLU;
    fpLU.setThreshold(1e-7);
    fpLU.compute(system_matrix);
    Eigen::VectorXd h_temp = fpLU.solve(u);
    cout << "Solution found" << endl;
    //cout <<std::setprecision(16)<< h_temp(783*3-1) << endl;

    /*Eigen::VectorXd h_temp = system_matrix.fullPivLu().solve(u);
    cout <<std::setprecision(16)<< h_temp(783*3-1) << endl;
    VectorXd test_vector = system_matrix*h_temp;*/

    std::ofstream file("H.txt");
    if (file.is_open())
    {
        file <<std::setprecision(16)<< h_temp << endl;
    }
    /*std::ofstream file4("test.txt");
    if (file4.is_open())
    {
        file4 << test_vector << endl;
    }*/
    //double relative_error = (system_matrix*h - u).norm() / u.norm();
    //std::cout << "The relative error is:\n" << relative_error << std::endl;
    //write the solutions to a text file
    
    /*std::ofstream file2("GS.txt");
    if (file2.is_open())
    {
        file2 << system_matrix << endl;
    }*/
    //cout << h << endl;
    return h;
    
}

MatrixX3d Greens_PWexpansion::solve4grid(std::string grid_path)
{
    cout << std::setprecision(16);
    MatrixXd grid = readMatrix(grid_path.c_str());
    int total_point_num = grid.rows();
    Matrix<double,Dynamic,Dynamic> us(total_point_num,3);
    cout << "Solver main loop starting" << endl;
    for (int i = 0; i < total_point_num; i++)
    {
        MatrixXd G_row(3,this->total_atom_num*3);
        Vector3d t_pos;
        t_pos(0) = grid(i,0);
        t_pos(1) = grid(i,1);
        t_pos(2) = grid(i,2);
        for (int atom_ind = 0; atom_ind < this->total_atom_num; atom_ind++)
        {
            Matrix3d G_partial = this->getGreens(atom_ind,t_pos(0),t_pos(1),t_pos(2));
            G_row.block(0,atom_ind*3,3,3) = G_partial;
        }
        Vector3d u_partial = G_row * this->h;
        us(i,0) = u_partial(0);
        us(i,1) = u_partial(1);
        us(i,2) = u_partial(2);
    }
    
    cout << "loop ended" << endl;
    //cout << us << endl;
    return us;
}

void Greens_PWexpansion::construct_eigenstrain_mat(){
    //this function constructs gradient of G tensor.
    //gradient components of G are augmented along rows
    cout << "beginning of construct_eigenstrain_mat" << endl;
    //get the displacements we are interested in
    VectorXd u(this->total_atom_num*3);
    for (int i = 0; i < this->total_atom_num; i++)
    {
        u(i*3) = this->disps(i,4);
        u(i*3+1) = this->disps(i,5);
        u(i*3+2) = this->disps(i,6);
    }
    
    std::ofstream file3("solve4.txt");
    if (file3.is_open())
    {
        file3 << u << endl;
    }

    //build the system matrix
    MatrixXd grad_matrix(this->total_atom_num*3*3, this->total_atom_num*3);
    cout << "system matrix loop begins" << endl;
    auto start = chrono::high_resolution_clock::now();
    double perturbation = 0.001;
    for (int grad_ind = 0; grad_ind < 3; grad_ind++)
    {
        MatrixXd system_matrix(this->total_atom_num*3, this->total_atom_num*3);

        for (int t_ind = 0; t_ind < this->total_atom_num; t_ind++)
        {
            Vector3d t_pos;
            t_pos(0) = this->disps(t_ind,1);
            t_pos(1) = this->disps(t_ind,2);
            t_pos(2) = this->disps(t_ind,3);
            for (int s_ind = 0; s_ind < this->total_atom_num; s_ind++)
            {
                t_pos(grad_ind) = this->disps(t_ind,grad_ind) + perturbation;
                Matrix3d G_partial_forward = this->getGreens(s_ind,t_pos(0), t_pos(1), t_pos(2));
                t_pos(grad_ind) = this->disps(t_ind,grad_ind) - 2 * perturbation;
                Matrix3d G_partial_backward = this->getGreens(s_ind,t_pos(0), t_pos(1), t_pos(2));
                Matrix3d G_partial = (G_partial_forward - G_partial_backward) / (2 * perturbation);
                system_matrix.block(t_ind*3,s_ind*3,3,3) = G_partial;
            }
            
        }
        cout << "Grad matrix:" << endl;
        grad_matrix.block(grad_ind * this->total_atom_num * 3, 0, this->total_atom_num*3, this->total_atom_num*3) = system_matrix;
    }
    
    
    
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "system matrix is created in " << duration.count() << " ms" << endl;
    
    //write the solutions to a text file
    
    std::ofstream file2("G_GRAD.txt");
    if (file2.is_open())
    {
        file2 << grad_matrix << endl;
    }
}

void Greens_PWexpansion::set_dislocation_per_vec(double v1, double v2, double v3){
    //set the periodicity vector along the dislocation line
    this->disloc_per_vec(0) = v1;
    this->disloc_per_vec(1) = v2;
    this->disloc_per_vec(2) = v3;
}

Matrix3d Greens_PWexpansion::getDislocGreens(int atom_index, double px, double py, double pz, int max_image_num){
    //this function gets the total Green's function for an atom in dislocation configuration
    //only finite number of images are considered along dislocation line (-max_image_num,+max_image_num)
    //cout << "index: " << atom_index <<endl;
    Eigen::Matrix3d G_partial = this->getGreens(atom_index,px,py,pz);
    //cout << "debug" << endl;
    if (max_image_num == 0){
        return G_partial;
    }
    
    for (int image_num = 1; image_num < max_image_num + 1; image_num++){
        //positive image
        Vector3d t_pos;
        t_pos(0) = px + image_num * this->disloc_per_vec(0);
        t_pos(1) = py + image_num * this->disloc_per_vec(1);
        t_pos(2) = pz + image_num * this->disloc_per_vec(2);
        G_partial = G_partial + this->getGreens(atom_index,t_pos(0),t_pos(1),t_pos(2));
        //negative image
        t_pos(0) = px - image_num * this->disloc_per_vec(0);
        t_pos(1) = py - image_num * this->disloc_per_vec(1);
        t_pos(2) = pz - image_num * this->disloc_per_vec(2);
        G_partial = G_partial + this->getGreens(atom_index,t_pos(0),t_pos(1),t_pos(2));
    }
    return G_partial;
}

VectorXd Greens_PWexpansion::solve_h_disloc(int max_image_num){
    //this function finds an h vector. However, when it writes it, it is just bogus. I could not solve for the life of me
    cout << "beginning of solve_h" << endl;
    //get the displacements we are interested in
    VectorXd u(this->total_atom_num*3);
    for (int i = 0; i < this->total_atom_num; i++)
    {
        u(i*3) = this->disps(i,4);
        u(i*3+1) = this->disps(i,5);
        u(i*3+2) = this->disps(i,6);
    }
    
    std::ofstream file3("solve4.txt");
    if (file3.is_open())
    {
        file3 << u << endl;
    }
    //cout << "original:" << endl << this->disps.block(0,4,this->total_atom_num,3) << endl;
    //cout << "reshaped:" << endl << this->disps.block(0,4,this->total_atom_num,3).reshaped<RowMajor>() << endl;
    //cout << u << endl;

    //build the system matrix
    MatrixXd system_matrix(this->total_atom_num*3, this->total_atom_num*3);
    cout << "system matrix loop begins" << endl;
    auto start = chrono::high_resolution_clock::now();
    for (int t_ind = 0; t_ind < this->total_atom_num; t_ind++)
    {
        Vector3d t_pos;
        t_pos(0) = this->disps(t_ind,1);
        t_pos(1) = this->disps(t_ind,2);
        t_pos(2) = this->disps(t_ind,3);
        cout << "t: " << t_ind << endl;
        //for (int s_ind = 0; s_ind < this->total_atom_num; s_ind++)
        //{
        //    //cout << "s: " << s_ind << endl;
        //    Matrix3d G_partial = this->getDislocGreens(s_ind,t_pos(0), t_pos(1), t_pos(2),max_image_num);
        //    system_matrix.block(t_ind*3,s_ind*3,3,3) = G_partial;
        //    
        //}
        this->async_par_for(0,this->total_atom_num,[&](unsigned s){
                          this->sys_mat_row_par(s,t_ind,t_pos(0), t_pos(1), t_pos(2),max_image_num,system_matrix);   //do something here with the index i
                        });
        
    }
    
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "system matrix is created in " << duration.count() << " ms" << endl;
    cout << "system matrix shape rows x cols: " << system_matrix.rows() << " x " << system_matrix.cols() << endl;
    cout << "u shape: " << u.rows() << " x " << u.cols() << endl;
    //solve the linear system for the fitting coefficients
    
    Eigen::FullPivLU<MatrixXd> fpLU;
    fpLU.setThreshold(1e-7);
    fpLU.compute(system_matrix);
    Eigen::VectorXd h_temp = fpLU.solve(u);
    cout << "Solution found" << endl;
    cout << "h shape: " << h_temp.rows() << " x " << h_temp.cols() << endl;
    //cout << "h: " << endl;
    //cout << h_temp << endl;
    //cout <<std::setprecision(16)<< h_temp(0) << endl;
    //cout <<std::setprecision(16)<< h_temp(1) << endl;
    //cout <<std::setprecision(16)<< h_temp(2) << endl;

    /*Eigen::VectorXd h_temp = system_matrix.fullPivLu().solve(u);
    cout <<std::setprecision(16)<< h_temp(783*3-1) << endl;*/
    //Eigen::VectorXd test_vector = system_matrix*h_temp;

    std::ofstream file("H.txt");
    if (file.is_open())
    {
        file <<std::setprecision(16)<< h_temp << endl;
    }
    /*std::ofstream file4("test.txt");
    if (file4.is_open())
    {
        file4 << test_vector << endl;
    }*/
    //double relative_error = (system_matrix*h_temp - u).norm() / u.norm();
    //std::cout << "The relative error is:\n" << relative_error << std::endl;
    //write the solutions to a text file
    
    /*std::ofstream file2("GS.txt");
    if (file2.is_open())
    {
        file2 << system_matrix << endl;
    }*/
    //cout << h << endl;
    return h;
}

MatrixX3d Greens_PWexpansion::solve4grid_disloc(std::string grid_path, int max_image_num){
    cout << std::setprecision(16);
    MatrixXd grid = readMatrix(grid_path.c_str());
    int total_point_num = grid.rows();
    Matrix<double,Dynamic,Dynamic> us(total_point_num,3);
    cout << "Solver main loop starting" << endl;
    for (int i = 0; i < total_point_num; i++)
    {
        MatrixXd G_row(3,this->total_atom_num*3);
        Vector3d t_pos;
        t_pos(0) = grid(i,0);
        t_pos(1) = grid(i,1);
        t_pos(2) = grid(i,2);
        for (int atom_ind = 0; atom_ind < this->total_atom_num; atom_ind++)
        {
            Matrix3d G_partial = this->getDislocGreens(atom_ind,t_pos(0),t_pos(1),t_pos(2),max_image_num);
            G_row.block(0,atom_ind*3,3,3) = G_partial;
        }
        Vector3d u_partial = G_row * this->h;
        us(i,0) = u_partial(0);
        us(i,1) = u_partial(1);
        us(i,2) = u_partial(2);
    }
    
    cout << "loop ended" << endl;
    //cout << us << endl;
    return us;
}

MatrixX3d Greens_PWexpansion::solve4grid_disloc_threaded(std::string grid_path, int max_image_num){
    cout << std::setprecision(16);
    cout << "grid path: " << grid_path << endl;
    MatrixXd grid = readMatrix(grid_path.c_str());
    int total_point_num = grid.rows();
    Matrix<double,Dynamic,Dynamic> us(total_point_num,3);
    cout << "Solver main loop starting" << endl;
    for (int i = 0; i < total_point_num; i++)
    {
        MatrixXd G_row(3,this->total_atom_num*3);
        Vector3d t_pos;
        t_pos(0) = grid(i,0);
        t_pos(1) = grid(i,1);
        t_pos(2) = grid(i,2);
        //for (int atom_ind = 0; atom_ind < this->total_atom_num; atom_ind++)
        //{
        //    Matrix3d G_partial = this->getDislocGreens(atom_ind,t_pos(0),t_pos(1),t_pos(2),max_image_num);
        //    G_row.block(0,atom_ind*3,3,3) = G_partial;
        //}
        this->async_par_for(0,this->total_atom_num,[&](unsigned s){
                          this->G_row_par(s,t_pos(0), t_pos(1), t_pos(2),max_image_num,G_row);   //do something here with the index i
                        });
        Vector3d u_partial = G_row * this->h;
        //cout << "u: " << u_partial << endl;
        us(i,0) = u_partial(0);
        us(i,1) = u_partial(1);
        us(i,2) = u_partial(2);
    }
    
    cout << "loop ended" << endl;
    //cout << us << endl;
    return us;
}

void Greens_PWexpansion::G_row_par(int s, double x, double y, double z, int max_image_num, Eigen::Ref<Eigen::MatrixXd> G_row){
    Matrix3d G_partial = this->getDislocGreens(s,x,y,z,max_image_num);
    //cout << "s: " << s << endl;
    //cout << "G: " << endl;
    //cout << G_partial << endl;
    //sys_mat(t*3+0,s*3+0) = G_partial(0,0);
    G_row.block(0,s*3,3,3) = G_partial;
}

void Greens_PWexpansion::sys_mat_row_par(int s, int t, double x, double y, double z, int max_image_num, Eigen::Ref<Eigen::MatrixXd> sys_mat){
    //cout << "debug sys mat assignment: " << s << endl;
    Matrix3d G_partial = this->getDislocGreens(s,x,y,z,max_image_num);
    //sys_mat(t*3+0,s*3+0) = G_partial(0,0);
    sys_mat.block(t*3,s*3,3,3) = G_partial;
}

void Greens_PWexpansion::async_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn){
    //internal loop
        bool par = true;
        auto int_fn = [&fn](unsigned int_start, unsigned seg_size){
            for (unsigned j = int_start; j < int_start+seg_size; j++){
                fn(j);
            }
        };

        //sequenced for
        if(!par){
            return int_fn(start, end);
        }

        //get number of threads
        unsigned nb_threads_hint = std::thread::hardware_concurrency();
        unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
        if(this->thread_num_set == false){
            //std::cout << "num threads: " << nb_threads << std::endl;
        }
        else{
            nb_threads = this->thread_num;
            //std::cout << "num threads: " << nb_threads << std::endl;
        }
        

        //calculate segments
        unsigned total_length = end - start;
        unsigned seg = total_length/nb_threads;
        unsigned last_seg = seg + total_length%nb_threads;
        
        //launch threads - parallel for
        auto fut_vec = std::vector<std::future<void>>();
        fut_vec.reserve(nb_threads);
        for(int k = 0; k < nb_threads-1; ++k){
            unsigned current_start = seg*k;
            fut_vec.emplace_back(std::async(int_fn, current_start, seg));
        }
        {
            unsigned current_start = seg*(nb_threads-1);
            fut_vec.emplace_back(std::async(std::launch::async, int_fn, current_start, last_seg));
        }
        for (auto& th : fut_vec){
            th.get();
        }
}

void Greens_PWexpansion::thread_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn){

        //internal loop
        bool par = true;
        auto int_fn = [&fn](unsigned int_start, unsigned seg_size){
            for (unsigned j = int_start; j < int_start+seg_size; j++){
                fn(j);
            }
        };

        //sequenced for
        if(!par){
            return int_fn(start, end);
        }

        //get number of threads
        unsigned nb_threads_hint = std::thread::hardware_concurrency();
        unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
        if(this->thread_num_set == false){
            //std::cout << "num threads: " << nb_threads << std::endl;
        }
        else{
            nb_threads = this->thread_num;
            //std::cout << "num threads: " << nb_threads << std::endl;
        }

        //calculate segments
        unsigned total_length = end - start;
        unsigned seg = total_length/nb_threads;
        unsigned last_seg = seg + total_length%nb_threads;

        //launch threads - parallel for
        auto threads_vec = std::vector<std::thread>();
        threads_vec.reserve(nb_threads);
        for(int k = 0; k < nb_threads-1; ++k){
            unsigned current_start = seg*k;
            threads_vec.emplace_back(std::thread(int_fn, current_start, seg));
        }
        {
            unsigned current_start = seg*(nb_threads-1);
            threads_vec.emplace_back(std::thread(int_fn, current_start, last_seg));
        }
        for (auto& th : threads_vec){
            th.join();
        }
    }

void Greens_PWexpansion::set_thread_num(int num){
        this->thread_num_set = true;
        this->thread_num = num;
    }