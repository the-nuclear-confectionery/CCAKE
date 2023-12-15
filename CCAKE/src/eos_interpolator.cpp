#include "eos_interpolator.h"

using namespace ccake;

EoS_Interpolator::EoS_Interpolator(fs::path path_to_eos_table)
{
  if (fs::exists(path_to_eos_table))
  {

    //Open hdf5 file for reading
    H5::H5File eos_file;
    try{
      eos_file = H5::H5File(path_to_eos_table.string(), H5F_ACC_RDONLY);
    } catch (H5::FileIException error){
      std::cout << "Error opening EoS table: " << error.getCDetailMsg() << std::endl;
      std::cout << "Exiting..." << std::endl;
      exit(1);
    }
    ///Retrieve the table information from the attributes
    H5::Attribute s_max_attr = eos_file.openAttribute("s_max");
    H5::Attribute ds_attr = eos_file.openAttribute("ds");
    H5::Attribute B_max_attr = eos_file.openAttribute("rhoB_max");
    H5::Attribute dB_attr = eos_file.openAttribute("drhoB");
    H5::Attribute Q_max_attr = eos_file.openAttribute("rhoQ_max");
    H5::Attribute dQ_attr = eos_file.openAttribute("drhoQ");
    H5::Attribute S_max_attr = eos_file.openAttribute("rhoS_max");
    H5::Attribute dS_attr = eos_file.openAttribute("drhoS");

    s_max_attr.read(H5::PredType::NATIVE_DOUBLE, &s_max);
    ds_attr.read(H5::PredType::NATIVE_DOUBLE, &ds);
    B_max_attr.read(H5::PredType::NATIVE_DOUBLE, &B_max);
    dB_attr.read(H5::PredType::NATIVE_DOUBLE, &dB);
    Q_max_attr.read(H5::PredType::NATIVE_DOUBLE, &Q_max);
    dQ_attr.read(H5::PredType::NATIVE_DOUBLE, &dQ);
    S_max_attr.read(H5::PredType::NATIVE_DOUBLE, &S_max);
    dS_attr.read(H5::PredType::NATIVE_DOUBLE, &dS);

    //Retrieve the axis information from the attributes
    H5::Attribute sAxis_attr = eos_file.openAttribute("s_Axis");
    H5::Attribute rhoBAxis_attr = eos_file.openAttribute("rhoB_Axis");
    H5::Attribute rhoQAxis_attr = eos_file.openAttribute("rhoQ_Axis");
    H5::Attribute rhoSAxis_attr = eos_file.openAttribute("rhoS_Axis");

    int s_Axis, rhoB_Axis, rhoQ_Axis, rhoS_Axis;
    sAxis_attr.read(H5::PredType::NATIVE_INT, &s_Axis);
    rhoBAxis_attr.read(H5::PredType::NATIVE_INT, &rhoB_Axis);
    rhoQAxis_attr.read(H5::PredType::NATIVE_INT, &rhoQ_Axis);
    rhoSAxis_attr.read(H5::PredType::NATIVE_INT, &rhoS_Axis);

    //Look at the temperature dataset to determine the shape of the table
    H5::DataSet T_dataset = eos_file.openDataSet("T");
    H5::DataSpace T_dataspace = T_dataset.getSpace();
    int T_rank = T_dataspace.getSimpleExtentNdims();
    hsize_t T_dims_out[T_rank];
    T_dataspace.getSimpleExtentDims(T_dims_out, NULL);
    Ns = T_dims_out[s_Axis];
    NB = T_dims_out[rhoB_Axis];
    NQ = T_dims_out[rhoQ_Axis];
    NS = T_dims_out[rhoS_Axis];
    //T_dataset.close();
    formatted_output::detail("s_max = " + std::to_string(s_max)
                            + ", ds = " + std::to_string(ds)
                            + ", Ns = " + std::to_string(Ns));
    formatted_output::detail("B_max = " + std::to_string(B_max)
                            + ", dB = " + std::to_string(dB)
                            + ", NB = " + std::to_string(NB));
    formatted_output::detail("Q_max = " + std::to_string(Q_max)
                            + ", dQ = " + std::to_string(dQ)
                            + ", NQ = " + std::to_string(NQ));
    formatted_output::detail("S_max = " + std::to_string(S_max)
                            + ", dS = " + std::to_string(dS)
                            + ", NS = " + std::to_string(NS));

    //Allocate memory for the thermodynamic variables in the device

    // T = eos_thermo_nonconst("T", Ns, NB, NQ, NS);
    // muB = eos_thermo_nonconst("muB", Ns, NB, NQ, NS);
    // muQ = eos_thermo_nonconst("muQ", Ns, NB, NQ, NS);
    // muS = eos_thermo_nonconst("muS", Ns, NB, NQ, NS);
    // e = eos_thermo_nonconst("e", Ns, NB, NQ, NS);
    // p = eos_thermo_nonconst("p", Ns, NB, NQ, NS);
    // cs2 = eos_thermo_nonconst("cs2", Ns, NB, NQ, NS);
    // dw_ds = eos_thermo_nonconst("dw_ds", Ns, NB, NQ, NS);
    // dw_dB = eos_thermo_nonconst("dw_dB", Ns, NB, NQ, NS);
    // dw_dQ = eos_thermo_nonconst("dw_dQ", Ns, NB, NQ, NS);
    // dw_dS = eos_thermo_nonconst("dw_dS", Ns, NB, NQ, NS);

    //Create a mirror view of the thermodynamic variables in the host
    // auto T_host = Kokkos::create_mirror_view(T);
    // auto muB_host = Kokkos::create_mirror_view(muB);
    // auto muQ_host = Kokkos::create_mirror_view(muQ);
    // auto muS_host = Kokkos::create_mirror_view(muS);
    // auto e_host = Kokkos::create_mirror_view(e);
    // auto p_host = Kokkos::create_mirror_view(p);
    // auto cs2_host = Kokkos::create_mirror_view(cs2);
    // auto dw_ds_host = Kokkos::create_mirror_view(dw_ds);
    // auto dw_dB_host = Kokkos::create_mirror_view(dw_dB);
    // auto dw_dQ_host = Kokkos::create_mirror_view(dw_dQ);
    // auto dw_dS_host = Kokkos::create_mirror_view(dw_dS);


    eos_vars[ccake::eos_variables::T] = eos_thermo_nonconst("T", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::muB] = eos_thermo_nonconst("muB", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::muQ] = eos_thermo_nonconst("muQ", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::muS] = eos_thermo_nonconst("muS", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::e] = eos_thermo_nonconst("e", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::p] = eos_thermo_nonconst("p", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::cs2] = eos_thermo_nonconst("cs2", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::dw_ds] = eos_thermo_nonconst("dw_ds", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::dw_dB] = eos_thermo_nonconst("dw_dB", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::dw_dQ] = eos_thermo_nonconst("dw_dQ", Ns, NB, NQ, NS);
    eos_vars[ccake::eos_variables::dw_dS] = eos_thermo_nonconst("dw_dS", Ns, NB, NQ, NS);

    auto T_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::T]);
    auto muB_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::muB]);
    auto muQ_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::muQ]);
    auto muS_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::muS]);
    auto e_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::e]);
    auto p_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::p]);
    auto cs2_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::cs2]);
    auto dw_ds_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::dw_ds]);
    auto dw_dB_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::dw_dB]);
    auto dw_dQ_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::dw_dQ]);
    auto dw_dS_host = Kokkos::create_mirror_view(eos_vars[ccake::eos_variables::dw_dS]);

    //Open datasets
    //H5::DataSet T_dataset = eos_file.openDataSet("T"); //Already opened above
    H5::DataSet muB_dataset = eos_file.openDataSet("muB");
    H5::DataSet muQ_dataset = eos_file.openDataSet("muQ");
    H5::DataSet muS_dataset = eos_file.openDataSet("muS");
    H5::DataSet e_dataset = eos_file.openDataSet("e");
    H5::DataSet p_dataset = eos_file.openDataSet("p");
    H5::DataSet cs2_dataset = eos_file.openDataSet("cs2");
    H5::DataSet dw_ds_dataset = eos_file.openDataSet("dwds");
    H5::DataSet dw_dB_dataset = eos_file.openDataSet("dwdB");
    H5::DataSet dw_dQ_dataset = eos_file.openDataSet("dwdQ");
    H5::DataSet dw_dS_dataset = eos_file.openDataSet("dwdS");

    //Create buffers for reading datasets
    std::vector<double> T_buffer(Ns*NB*NQ*NS);
    std::vector<double> muB_buffer(Ns*NB*NQ*NS);
    std::vector<double> muQ_buffer(Ns*NB*NQ*NS);
    std::vector<double> muS_buffer(Ns*NB*NQ*NS);
    std::vector<double> e_buffer(Ns*NB*NQ*NS);
    std::vector<double> p_buffer(Ns*NB*NQ*NS);
    std::vector<double> cs2_buffer(Ns*NB*NQ*NS);
    std::vector<double> dw_ds_buffer(Ns*NB*NQ*NS);
    std::vector<double> dw_dB_buffer(Ns*NB*NQ*NS);
    std::vector<double> dw_dQ_buffer(Ns*NB*NQ*NS);
    std::vector<double> dw_dS_buffer(Ns*NB*NQ*NS);

    //Read datasets
    T_dataset.read(T_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    muB_dataset.read(muB_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    muQ_dataset.read(muQ_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    muS_dataset.read(muS_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    e_dataset.read(e_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    p_dataset.read(p_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    cs2_dataset.read(cs2_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    dw_ds_dataset.read(dw_ds_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    dw_dB_dataset.read(dw_dB_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    dw_dQ_dataset.read(dw_dQ_buffer.data(), H5::PredType::NATIVE_DOUBLE);
    dw_dS_dataset.read(dw_dS_buffer.data(), H5::PredType::NATIVE_DOUBLE);

    //Transfer buffers to the Kokkos host views
    assert(s_Axis == 0 && rhoB_Axis == 1 && rhoS_Axis == 2 && rhoQ_Axis == 3);
    for(int i = 0; i < Ns*NB*NQ*NS ; i++)
    {
      int is = i/(NB*NQ*NS);
      int iB = (i - is*NB*NQ*NS)/(NQ*NS);
      int iQ = (i - is*NB*NQ*NS - iB*NQ*NS)/NS;
      int iS = i - is*NB*NQ*NS - iB*NQ*NS - iQ*NS;

      T_host(is, iB, iS, iQ) = T_buffer[i];
      muB_host(is, iB, iS, iQ) = muB_buffer[i];
      muQ_host(is, iB, iS, iQ) = muQ_buffer[i];
      muS_host(is, iB, iS, iQ) = muS_buffer[i];
      e_host(is, iB, iS, iQ) = e_buffer[i];
      p_host(is, iB, iS, iQ) = p_buffer[i];
      cs2_host(is, iB, iS, iQ) = cs2_buffer[i];
      dw_ds_host(is, iB, iS, iQ) = dw_ds_buffer[i];
      dw_dB_host(is, iB, iS, iQ) = dw_dB_buffer[i];
      dw_dQ_host(is, iB, iS, iQ) = dw_dQ_buffer[i];
      dw_dS_host(is, iB, iS, iQ) = dw_dS_buffer[i];
    }

    // Kokkos::deep_copy(T, T_host);
    // Kokkos::deep_copy(muB, muB_host);
    // Kokkos::deep_copy(muQ, muQ_host);
    // Kokkos::deep_copy(muS, muS_host);
    // Kokkos::deep_copy(e, e_host);
    // Kokkos::deep_copy(p, p_host);
    // Kokkos::deep_copy(cs2, cs2_host);
    // Kokkos::deep_copy(dw_ds, dw_ds_host);
    // Kokkos::deep_copy(dw_dB, dw_dB_host);
    // Kokkos::deep_copy(dw_dQ, dw_dQ_host);
    // Kokkos::deep_copy(dw_dS, dw_dS_host);
//
    //Lastly, init the const views
    // eos_vars[ccake::eos_variables::T] = new eos_thermo(T);
    // eos_vars[ccake::eos_variables::muB] = new eos_thermo(muB);
    // eos_vars[ccake::eos_variables::muQ] = new eos_thermo(muQ);
    // eos_vars[ccake::eos_variables::muS] = new eos_thermo(muS);
    // eos_vars[ccake::eos_variables::e] = new eos_thermo(e);
    // eos_vars[ccake::eos_variables::p] = new eos_thermo(p);
    // eos_vars[ccake::eos_variables::cs2] = new eos_thermo(cs2);
    // eos_vars[ccake::eos_variables::dw_ds] = new eos_thermo(dw_ds);
    // eos_vars[ccake::eos_variables::dw_dB] = new eos_thermo(dw_dB);
    // eos_vars[ccake::eos_variables::dw_dQ] = new eos_thermo(dw_dQ);
    // eos_vars[ccake::eos_variables::dw_dS] = new eos_thermo(dw_dS);

    Kokkos::deep_copy(eos_vars[ccake::eos_variables::T], T_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::muB], muB_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::muQ], muQ_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::muS], muS_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::e], e_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::p], p_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::cs2], cs2_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::dw_ds], dw_ds_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::dw_dB], dw_dB_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::dw_dQ], dw_dQ_host);
    Kokkos::deep_copy(eos_vars[ccake::eos_variables::dw_dS], dw_dS_host);

    formatted_output::detail("EoS table loaded from " + path_to_eos_table.string());
    eos_file.close();
  }
  else
  {
    std::cout << "EoS table not found at " << path_to_eos_table << std::endl;
    std::cout << "Exiting..." << std::endl;
    exit(1);
  }
}

void EoS_Interpolator::fill_thermodynamics(Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH> &particles, const double t){
  CREATE_VIEW(device_,particles)
  auto interpolate = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
    //Compute the entropy and charge densities in the particles' rest frame
    double s = device_thermo.access(is, ia, ccake::thermo_info::s);
    double rhoB = Kokkos::fabs(device_thermo.access(is, ia, ccake::thermo_info::rhoB));
    double rhoS = Kokkos::fabs(device_thermo.access(is, ia, ccake::thermo_info::rhoS));
    double rhoQ = Kokkos::fabs(device_thermo.access(is, ia, ccake::thermo_info::rhoQ));

    //Find the indices of the nearest neighbors in the table
    int idx[4];
    idx[0] = (int) Kokkos::floor(s/ds);
    idx[1] = (int) Kokkos::floor(rhoB/dB);
    idx[2] = (int) Kokkos::floor(rhoS/dS);
    idx[3] = (int) Kokkos::floor(rhoQ/dQ);

    double pos[4];
    pos[0] = s;
    pos[1] = rhoB;
    pos[2] = rhoS;
    pos[3] = rhoQ;

    device_thermo.access(is, ia, ccake::thermo_info::T)    = interpolate4D(idx, pos, ccake::eos_variables::T);
    device_thermo.access(is, ia, ccake::thermo_info::muB)  = interpolate4D(idx, pos, ccake::eos_variables::muB);
    device_thermo.access(is, ia, ccake::thermo_info::muQ)  = interpolate4D(idx, pos, ccake::eos_variables::muQ);
    device_thermo.access(is, ia, ccake::thermo_info::muS)  = interpolate4D(idx, pos, ccake::eos_variables::muS);
    device_thermo.access(is, ia, ccake::thermo_info::e)    = interpolate4D(idx, pos, ccake::eos_variables::e);
    device_thermo.access(is, ia, ccake::thermo_info::p)    = interpolate4D(idx, pos, ccake::eos_variables::p);
    device_thermo.access(is, ia, ccake::thermo_info::cs2)  = interpolate4D(idx, pos, ccake::eos_variables::cs2);
    device_thermo.access(is, ia, ccake::thermo_info::dwds) = interpolate4D(idx, pos, ccake::eos_variables::dw_ds);
    device_thermo.access(is, ia, ccake::thermo_info::dwdB) = interpolate4D(idx, pos, ccake::eos_variables::dw_dB);
    device_thermo.access(is, ia, ccake::thermo_info::dwdS) = interpolate4D(idx, pos, ccake::eos_variables::dw_dS);
    device_thermo.access(is, ia, ccake::thermo_info::dwdQ) = interpolate4D(idx, pos, ccake::eos_variables::dw_dQ);

    device_thermo.access(is, ia, ccake::thermo_info::w) = device_thermo.access(is, ia, ccake::thermo_info::e)
                                                        + device_thermo.access(is, ia, ccake::thermo_info::p);
    device_thermo.access(is, ia, ccake::thermo_info::A) = device_thermo.access(is, ia, ccake::thermo_info::w)
                                                         -s*device_thermo.access(is, ia, ccake::thermo_info::dwds);
  };
  Cabana::SimdPolicy<VECTOR_LENGTH, ExecutionSpace> simd_policy(0, particles.size());
  Cabana::simd_parallel_for(simd_policy, interpolate, "interpolate");
  Kokkos::fence();

} // namespace ccake

/// @brief Interpolates the EoS table
/// @details For the algorithm used, see https://jmlr.org/papers/volume17/15-243/15-243.pdf
/// section 3.1. https://dl.acm.org/doi/fullHtml/10.1145/3423184 discusses possible
/// improvements to the algorithm, which could be implemented in the future.
KOKKOS_FUNCTION
double EoS_Interpolator::interpolate4D(int idx [], double pos[], int ivar ) const{

  //Check if the point is outside the grid
  if(idx[0] < 0 || idx[0] >= Ns-1 ||
     idx[1] < 0 || idx[1] >= NB-1 ||
     idx[2] < 0 || idx[2] >= NS-1 ||
     idx[3] < 0 || idx[3] >= NQ-1)
  {
    return NAN; //I cannot print a warning here because this function is called from a device kernel.
                //Instead, NAN will ensure that the simulation crashes if outside the table.
  }
  ///Creates a vector x which lies in the unit hypercube
  double x[4];
  x[0] = (pos[0] - idx[0]*ds)/ds;
  x[1] = (pos[1] - idx[1]*dB)/dB;
  x[2] = (pos[2] - idx[2]*dS)/dS;
  x[3] = (pos[3] - idx[3]*dQ)/dQ;

  //Asserts that we are in the unit hypercube
  #ifdef DEBUG
  //for(int i = 0; i < 4; i++) assert(x[i] >= 0 && x[i] <= 1);
  #endif
  auto bit = [](int i, int k){return (k >> i) & 1;};
  double phi[16]; //A function phi for each vertex of the hypercube
  double theta[16]; //The value to be interpolates in each point of the hypercube
  for(int k = 0; k < 16; k++)
  {
    phi[k] = 1;
    for(int i = 0; i < 4; i++) phi[k] *= bit(i, k) ? x[i] : (1 - x[i]);
  }

  for(int k = 0; k < 16; ++k)
  {
    int idx_k[4];
    for(int i = 0; i < 4; ++i)
      idx_k[i] = idx[i] + bit(i, k);
    theta[k] = eos_vars[ivar](idx_k[0], idx_k[1], idx_k[2], idx_k[3]);
    //cout << idx_k[0] << " " << idx_k[1] << " " << idx_k[2] << " " << idx_k[3] << " " << phi[k] << " "  << theta[k] << endl;
  }
  double f = 0;

  for(int k = 0; k < 16; k++) f += phi[k]*theta[k];

  return f;
}
