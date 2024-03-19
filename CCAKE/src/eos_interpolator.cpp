///@file eos_interpolator
///@author Willian M. Serenone
///@date 2023-12-15
///@brief Implementation of the EoS interpolator class
#include "eos_interpolator.h"

using namespace ccake;

EoS_Interpolator::EoS_Interpolator(fs::path path_to_eos_table)
{
  if (fs::exists(path_to_eos_table))
  {

    //Open hdf5 file for reading
    try{
      eos_file = H5::H5File(path_to_eos_table.string(), H5F_ACC_RDONLY);
    } catch (H5::FileIException error){
      std::cout << "Error opening EoS table: " << error.getCDetailMsg() << std::endl;
      std::cout << "Exiting..." << std::endl;
      exit(1);
    }
    //Loop over tables, loading them into memory
    for (int is=-3; is<2; ++is)
    for (int iB=-3; iB<1; ++iB)
    for (int iS=-3; iS<1; ++iS)
    for (int iQ=-3; iQ<1; ++iQ)
      load_table({is, iB, iS, iQ});

    formatted_output::update("Loaded EoS to device memory");
    eos_file.close();
  }
  else
  {
    std::cout << "EoS table not found at " << path_to_eos_table << std::endl;
    std::cout << "Exiting..." << std::endl;
    exit(1);
  }
}

/// @brief Helper to load the header of one table of the EoS to memmory.
/// @details The table is identified by a list of exponents that indicates
/// the order of magnitude of the independent variables (s, rhoB, rhoS, rhoQ).
/// The list **must be** in the order (s, rhoB, rhoS, rhoQ) and have size 4.
/// @param exponents list of orders of magnitude of the independent variables.
void EoS_Interpolator::load_header(std::initializer_list<int> exponents){
  assert(exponents.size() == 4);
  std::string table_name = "table";
  for (auto exponent : exponents)
  {
    table_name += "_" + std::to_string(exponent);
  }

  //Get the attributes of the group "table_name"
  H5::Group group = eos_file.openGroup(table_name);
  int ts = *exponents.begin() + 3;
  int tB = *(exponents.begin() + 1) + 3;
  int tS = *(exponents.begin() + 2) + 3;
  int tQ = *(exponents.begin() + 3) + 3;

  H5::Attribute s_min_attr = group.openAttribute("s_min");
  H5::Attribute s_max_attr = group.openAttribute("s_max");
  H5::Attribute ds_attr = group.openAttribute("ds");
  H5::Attribute B_min_attr = group.openAttribute("rhoB_min");
  H5::Attribute B_max_attr = group.openAttribute("rhoB_max");
  H5::Attribute dB_attr = group.openAttribute("drhoB");
  H5::Attribute S_min_attr = group.openAttribute("rhoS_min");
  H5::Attribute S_max_attr = group.openAttribute("rhoS_max");
  H5::Attribute dS_attr = group.openAttribute("drhoS");
  H5::Attribute Q_min_attr = group.openAttribute("rhoQ_min");
  H5::Attribute Q_max_attr = group.openAttribute("rhoQ_max");
  H5::Attribute dQ_attr = group.openAttribute("drhoQ");

  s_min_attr.read(H5::PredType::NATIVE_DOUBLE, &s_attr[ts][tB][tS][tQ].min);
  s_max_attr.read(H5::PredType::NATIVE_DOUBLE, &s_attr[ts][tB][tS][tQ].max);
  ds_attr.read(H5::PredType::NATIVE_DOUBLE, &s_attr[ts][tB][tS][tQ].step);
  B_min_attr.read(H5::PredType::NATIVE_DOUBLE, &B_attr[ts][tB][tS][tQ].min);
  B_max_attr.read(H5::PredType::NATIVE_DOUBLE, &B_attr[ts][tB][tS][tQ].max);
  dB_attr.read(H5::PredType::NATIVE_DOUBLE, &B_attr[ts][tB][tS][tQ].step);
  S_min_attr.read(H5::PredType::NATIVE_DOUBLE, &S_attr[ts][tB][tS][tQ].min);
  S_max_attr.read(H5::PredType::NATIVE_DOUBLE, &S_attr[ts][tB][tS][tQ].max);
  dS_attr.read(H5::PredType::NATIVE_DOUBLE, &S_attr[ts][tB][tS][tQ].step);
  Q_min_attr.read(H5::PredType::NATIVE_DOUBLE, &Q_attr[ts][tB][tS][tQ].min);
  Q_max_attr.read(H5::PredType::NATIVE_DOUBLE, &Q_attr[ts][tB][tS][tQ].max);
  dQ_attr.read(H5::PredType::NATIVE_DOUBLE, &Q_attr[ts][tB][tS][tQ].step);

  //Retrieve the axis information from the attributes
  H5::Attribute sAxis_attr = group.openAttribute("s_Axis");
  H5::Attribute rhoBAxis_attr = group.openAttribute("rhoB_Axis");
  H5::Attribute rhoSAxis_attr = group.openAttribute("rhoS_Axis");
  H5::Attribute rhoQAxis_attr = group.openAttribute("rhoQ_Axis");

  auto get_N = [](double min, double max, double step){
    return (int) std::round(max/step);
  };

  //Compute the number of points in each direction
  s_attr[ts][tB][tS][tQ].N = get_N(s_attr[ts][tB][tS][tQ].min,
                                     s_attr[ts][tB][tS][tQ].max,
                                     s_attr[ts][tB][tS][tQ].step);
  B_attr[ts][tB][tS][tQ].N = get_N(B_attr[ts][tB][tS][tQ].min,
                                      B_attr[ts][tB][tS][tQ].max,
                                      B_attr[ts][tB][tS][tQ].step);
  S_attr[ts][tB][tS][tQ].N = get_N(S_attr[ts][tB][tS][tQ].min,
                                      S_attr[ts][tB][tS][tQ].max,
                                      S_attr[ts][tB][tS][tQ].step);
  Q_attr[ts][tB][tS][tQ].N = get_N(Q_attr[ts][tB][tS][tQ].min,
                                      Q_attr[ts][tB][tS][tQ].max,
                                      Q_attr[ts][tB][tS][tQ].step);

  //Asserts that the axis are in the correct order
  int s_Axis, rhoB_Axis, rhoQ_Axis, rhoS_Axis;
  sAxis_attr.read(H5::PredType::NATIVE_INT, &s_Axis);
  rhoBAxis_attr.read(H5::PredType::NATIVE_INT, &rhoB_Axis);
  rhoSAxis_attr.read(H5::PredType::NATIVE_INT, &rhoS_Axis);
  rhoQAxis_attr.read(H5::PredType::NATIVE_INT, &rhoQ_Axis);

  assert(s_Axis == 0 && rhoB_Axis == 1 && rhoS_Axis == 2 && rhoQ_Axis == 3);
  //Close group
  group.close();

}

/// @brief Load one table of the EoS to memmory.
/// @details The table is identified by a list of exponents that indicates
/// the order of magnitude of the independent variables (s, rhoB, rhoS, rhoQ).
/// The list **must be** in the order (s, rhoB, rhoS, rhoQ) and have size 4.
/// @param exponents list of orders of magnitude of the independent variables.
void EoS_Interpolator::load_table(std::initializer_list<int> exponents)
{
  assert(exponents.size() == 4);
  std::string table_name = "table";
  for (auto exponent : exponents)
  {
    table_name += "_" + std::to_string(exponent);
  }

  int ts = *exponents.begin() + 3;
  int tB = *(exponents.begin() + 1) + 3;
  int tS = *(exponents.begin() + 2) + 3;
  int tQ = *(exponents.begin() + 3) + 3;

  load_header(exponents);

  int Ns = s_attr[ts][tB][tS][tQ].N;
  int NB = B_attr[ts][tB][tS][tQ].N;
  int NS = S_attr[ts][tB][tS][tQ].N;
  int NQ = Q_attr[ts][tB][tS][tQ].N;

  //Open datasets
  H5::DataSet T_dataset = eos_file.openDataSet(table_name + "/T");
  H5::DataSet muB_dataset = eos_file.openDataSet(table_name + "/muB");
  H5::DataSet muS_dataset = eos_file.openDataSet(table_name + "/muS");
  H5::DataSet muQ_dataset = eos_file.openDataSet(table_name + "/muQ");
  H5::DataSet e_dataset = eos_file.openDataSet(table_name + "/e");
  H5::DataSet p_dataset = eos_file.openDataSet(table_name + "/p");
  H5::DataSet cs2_dataset = eos_file.openDataSet(table_name + "/cs2");
  H5::DataSet dw_ds_dataset = eos_file.openDataSet(table_name + "/dwds");
  H5::DataSet dw_dB_dataset = eos_file.openDataSet(table_name + "/dwdB");
  H5::DataSet dw_dS_dataset = eos_file.openDataSet(table_name + "/dwdS");
  H5::DataSet dw_dQ_dataset = eos_file.openDataSet(table_name + "/dwdQ");

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
  std::vector<double> dw_dS_buffer(Ns*NB*NQ*NS);
  std::vector<double> dw_dQ_buffer(Ns*NB*NQ*NS);

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
  dw_dS_dataset.read(dw_dS_buffer.data(), H5::PredType::NATIVE_DOUBLE);
  dw_dQ_dataset.read(dw_dQ_buffer.data(), H5::PredType::NATIVE_DOUBLE);

  //Allocate device memory
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::T] = Kokkos::View<double****,DeviceType>("T", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muB] = Kokkos::View<double****,DeviceType>("muB", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muS] = Kokkos::View<double****,DeviceType>("muS", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muQ] = Kokkos::View<double****,DeviceType>("muQ", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::e] = Kokkos::View<double****,DeviceType>("e", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::p] = Kokkos::View<double****,DeviceType>("p", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::cs2] = Kokkos::View<double****,DeviceType>("cs2", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_ds] = Kokkos::View<double****,DeviceType>("dw_ds", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dB] = Kokkos::View<double****,DeviceType>("dw_dB", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dS] = Kokkos::View<double****,DeviceType>("dw_dS", Ns, NB, NS, NQ);
  eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dQ] = Kokkos::View<double****,DeviceType>("dw_dQ", Ns, NB, NS, NQ);

  //Create kokkos host views
  auto T_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::T]);
  auto muB_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muB]);
  auto muS_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muS]);
  auto muQ_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muQ]);
  auto e_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::e]);
  auto p_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::p]);
  auto cs2_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::cs2]);
  auto dw_ds_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_ds]);
  auto dw_dB_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dB]);
  auto dw_dS_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dS]);
  auto dw_dQ_host = Kokkos::create_mirror_view(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dQ]);

  //Transfer buffers to the Kokkos host views
  for(int i = 0; i < Ns*NB*NQ*NS ; i++)
  {
    int is = i/(NB*NQ*NS);
    int iB = (i - is*NB*NQ*NS)/(NQ*NS);
    int iS = (i - is*NB*NS*NQ - iB*NS*NQ)/NQ;
    int iQ = i - is*NB*NS*NQ - iB*NS*NQ - iS*NQ;

    T_host(is, iB, iS, iQ) = T_buffer[i];
    muB_host(is, iB, iS, iQ) = muB_buffer[i];
    muQ_host(is, iB, iS, iQ) = muQ_buffer[i];
    muS_host(is, iB, iS, iQ) = muS_buffer[i];
    e_host(is, iB, iS, iQ) = e_buffer[i];
    p_host(is, iB, iS, iQ) = p_buffer[i];
    cs2_host(is, iB, iS, iQ) = cs2_buffer[i];
    dw_ds_host(is, iB, iS, iQ) = dw_ds_buffer[i];
    dw_dB_host(is, iB, iS, iQ) = dw_dB_buffer[i];
    dw_dS_host(is, iB, iS, iQ) = dw_dS_buffer[i];
    dw_dQ_host(is, iB, iS, iQ) = dw_dQ_buffer[i];
  }

    //Copy to device
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::T], T_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muB], muB_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muQ], muQ_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::muS], muS_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::e], e_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::p], p_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::cs2], cs2_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_ds], dw_ds_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dB], dw_dB_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dQ], dw_dQ_host);
    Kokkos::deep_copy(eos_vars[ts][tB][tS][tQ][ccake::eos_variables::dw_dS], dw_dS_host);

}

void EoS_Interpolator::fill_thermodynamics(Cabana::AoSoA<CabanaParticle,
                                                        DeviceType,
                                                        VECTOR_LENGTH> &particles,
                                                        const double t){
  CREATE_VIEW(device_,particles)
  auto interpolate = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
    //Compute the entropy and charge densities in the particles' rest frame
    double s = Kokkos::max(.001,device_thermo.access(is, ia, ccake::thermo_info::s));
    double rhoB = Kokkos::max(.001,Kokkos::fabs(device_thermo.access(is, ia, ccake::thermo_info::rhoB)));
    double rhoS = Kokkos::max(.001,Kokkos::fabs(device_thermo.access(is, ia, ccake::thermo_info::rhoS)));
    double rhoQ = Kokkos::max(.001,Kokkos::fabs(device_thermo.access(is, ia, ccake::thermo_info::rhoQ)));

    //Needs to find which table to use by finding the exponents of the independent variables
    auto get_exponent = [](double x){
      int ie = (int) Kokkos::floor(Kokkos::log10(Kokkos::fabs(x)) + 3);
      return Kokkos::max(ie, 0);
    };
    int ts = Kokkos::min(get_exponent(s),4);
    int tB = Kokkos::min(get_exponent(rhoB),3);
    int tS = Kokkos::min(get_exponent(rhoS),3);
    int tQ = Kokkos::min(get_exponent(rhoQ),3);

    int tble[4] = {ts, tB, tS, tQ};

    double s_min = s_attr[ts][tB][tS][tQ].min;
    double ds = s_attr[ts][tB][tS][tQ].step;
    double B_min = B_attr[ts][tB][tS][tQ].min;
    double dB = B_attr[ts][tB][tS][tQ].step;
    double S_min = S_attr[ts][tB][tS][tQ].min;
    double dS = S_attr[ts][tB][tS][tQ].step;
    double Q_min = Q_attr[ts][tB][tS][tQ].min;
    double dQ = Q_attr[ts][tB][tS][tQ].step;

    //Find the indices of the nearest neighbors in the table
    int idx[4];
    idx[0] = (int) Kokkos::floor((s-s_min)/ds);
    idx[1] = (int) Kokkos::floor((rhoB-B_min)/dB);
    idx[2] = (int) Kokkos::floor((rhoS-S_min)/dS);
    idx[3] = (int) Kokkos::floor((rhoQ-Q_min)/dQ);

    double pos[4];
    pos[0] = s;
    pos[1] = rhoB;
    pos[2] = rhoS;
    pos[3] = rhoQ;

    device_thermo.access(is, ia, ccake::thermo_info::T)    = interpolate4D(idx, pos, tble, ccake::eos_variables::T);
    device_thermo.access(is, ia, ccake::thermo_info::muB)  = interpolate4D(idx, pos, tble, ccake::eos_variables::muB);
    device_thermo.access(is, ia, ccake::thermo_info::muS)  = interpolate4D(idx, pos, tble, ccake::eos_variables::muS);
    device_thermo.access(is, ia, ccake::thermo_info::muQ)  = interpolate4D(idx, pos, tble, ccake::eos_variables::muQ);
    device_thermo.access(is, ia, ccake::thermo_info::e)    = interpolate4D(idx, pos, tble, ccake::eos_variables::e);
    device_thermo.access(is, ia, ccake::thermo_info::p)    = interpolate4D(idx, pos, tble, ccake::eos_variables::p);
    device_thermo.access(is, ia, ccake::thermo_info::cs2)  = interpolate4D(idx, pos, tble, ccake::eos_variables::cs2);
    device_thermo.access(is, ia, ccake::thermo_info::dwds) = interpolate4D(idx, pos, tble, ccake::eos_variables::dw_ds);
    device_thermo.access(is, ia, ccake::thermo_info::dwdB) = interpolate4D(idx, pos, tble, ccake::eos_variables::dw_dB);
    device_thermo.access(is, ia, ccake::thermo_info::dwdS) = interpolate4D(idx, pos, tble, ccake::eos_variables::dw_dS);
    device_thermo.access(is, ia, ccake::thermo_info::dwdQ) = interpolate4D(idx, pos, tble, ccake::eos_variables::dw_dQ);

    device_thermo.access(is, ia, ccake::thermo_info::w) = device_thermo.access(is, ia, ccake::thermo_info::e)
                                                        + device_thermo.access(is, ia, ccake::thermo_info::p);
    device_thermo.access(is, ia, ccake::thermo_info::A) = device_thermo.access(is, ia, ccake::thermo_info::w)
                                                         -s*device_thermo.access(is, ia, ccake::thermo_info::dwds);
  };
  Cabana::SimdPolicy<VECTOR_LENGTH, ExecutionSpace> simd_policy(0, particles.size());
  Cabana::simd_parallel_for(simd_policy, interpolate, "interpolate");
  Kokkos::fence();

} // namespace ccake

/// @brief Linearly interpolate one the thermodynamic quantities.
/// @details It needs the desired position and the lower corner of the hypercube.
/// Based on this, it creates a unit hypercube and evaluates the chosen variable
/// at its 16 vertices. It then
/// @details
// KOKKOS_FUNCTION
// double EoS_Interpolator::interpolate4D_slow(int idx [], double pos[], int ivar ) const{

//    int N[4] = {Ns, NB, NS, NQ};
//    //Avoid out of lower bounds
//    for(int i = 0; i < 4; i++)
//     idx[i] = Kokkos::max(idx[i], 0);
//   //Avoid upper bounds
//   for(int i = 0; i < 4; i++)
//     idx[i] = Kokkos::min(idx[i], N[i]-2);

//    ///Creates a vector x which lies in the unit hypercube
//    double x[4];
//    x[0] = (pos[0] - s_min - idx[0]*ds)/ds;
//    x[1] = (pos[1] - B_min - idx[1]*dB)/dB;
//    x[2] = (pos[2] - S_min - idx[2]*dS)/dS;
//    x[3] = (pos[3] - Q_min - idx[3]*dQ)/dQ;

//   //Create vertex points
//   std::array<double[4], 16> vertex;
//   std::array<double, 16> values;
//   for(int it = 0; it < 2; it++)
//   for(int ix = 0; ix < 2; ix++)
//   for(int iy = 0; iy < 2; iy++)
//   for(int iz = 0; iz < 2; iz++)
//   {
//     int k = it*8 + ix*4 + iy*2 + iz;
//     vertex[k][0] = it; vertex[k][1] = ix;
//     vertex[k][2] = iy; vertex[k][3] = iz;
//     values[k] = eos_vars[ivar](idx[0] + it,
//                                idx[1] + ix,
//                                idx[2] + iy,
//                                idx[3] + iz);
//   }

//   //Interpolate along the z axis
//   double interp_vals[8];
//   for(int it = 0; it < 2; it++)
//   for(int ix = 0; ix < 2; ix++)
//   for(int iy = 0; iy < 2; iy++)
//   {
//     int l = it*4 + ix*2 + iy;
//     int k1 = it*8 + ix*4 + iy*2;
//     int k2 = k1+1;
//     interp_vals[l] = interpolate1D(x[3], vertex[k1][3], vertex[k2][3],
//         	                          values[k1], values[k2]);
//   }

//   //Interpolate along the y axis
//   for(int it = 0; it < 2; it++)
//   for(int ix = 0; ix < 2; ix++)
//   {
//     int l = it*2 + ix;
//     int k1 = it*8 + ix*4;
//     int k2 = k1+2;
//     interp_vals[l] = interpolate1D(x[2], vertex[k1][2], vertex[k2][2],
//         	                          interp_vals[2*l], interp_vals[2*l+1]);
//   }

//   //Interpolate along the x axis
//   for(int it = 0; it < 2; it++)
//   {
//     int k1 = it*8;
//     int k2 = k1+4;
//     interp_vals[it] = interpolate1D(x[1], vertex[k1][1], vertex[k2][1],
//         	                          interp_vals[2*it], interp_vals[2*it+1]);
//   }

//   return interpolate1D(x[0], vertex[0][0], vertex[0][8],
//         	                          interp_vals[0], interp_vals[1]);

// }

// KOKKOS_FUNCTION
// double EoS_Interpolator::interpolate1D(double x, double x0, double x1,
//                                                 double y0, double y1) const
// {
//    return y0 + (y1 - y0)/(x1 - x0)*(x - x0);
// }

/// @brief Interpolates the EoS table
/// @details For the algorithm used, see https://jmlr.org/papers/volume17/15-243/15-243.pdf
/// section 3.1. https://dl.acm.org/doi/fullHtml/10.1145/3423184 discusses possible
/// improvements to the algorithm, which could be implemented in the future.
KOKKOS_FUNCTION
double EoS_Interpolator::interpolate4D(int idx [], double pos[],
                                       int tble[], int ivar ) const{

  double s_min = s_attr[tble[0]][tble[1]][tble[2]][tble[3]].min;
  double ds    = s_attr[tble[0]][tble[1]][tble[2]][tble[3]].step;
  double B_min = B_attr[tble[0]][tble[1]][tble[2]][tble[3]].min;
  double dB    = B_attr[tble[0]][tble[1]][tble[2]][tble[3]].step;
  double S_min = S_attr[tble[0]][tble[1]][tble[2]][tble[3]].min;
  double dS    = S_attr[tble[0]][tble[1]][tble[2]][tble[3]].step;
  double Q_min = Q_attr[tble[0]][tble[1]][tble[2]][tble[3]].min;
  double dQ    = Q_attr[tble[0]][tble[1]][tble[2]][tble[3]].step;

  ///Creates a vector x which lies in the unit hypercube
  double x[4];
  x[0] = (pos[0] - s_min - idx[0]*ds)/ds;
  x[1] = (pos[1] - B_min - idx[1]*dB)/dB;
  x[2] = (pos[2] - S_min - idx[2]*dS)/dS;
  x[3] = (pos[3] - Q_min - idx[3]*dQ)/dQ;

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
    theta[k] = eos_vars[tble[0]][tble[1]][tble[2]][tble[3]][ivar](idx_k[0], idx_k[1], idx_k[2], idx_k[3]);
  }
  double f = 0;

  for(int k = 0; k < 16; k++) f += phi[k]*theta[k];

  return f;
}
