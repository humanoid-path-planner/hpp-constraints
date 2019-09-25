#include <boost/assign/list_of.hpp>

#include <hpp/pinocchio/util.hh>
#include <hpp/constraints/symbolic-function.hh>
#include <hpp/constraints/solver/hierarchical-iterative.hh>
#include <sstream>

using namespace hpp::constraints;

// This unit test was used to write function hpp::fcl::details::GJK::projectTetrahedraOrigin.

#define DECLARE_EXPRESSION(Class) \
  typedef typename Traits <Class>::Ptr_t Ptr_t;

#define HPP_CONSTRAINTS_CB_CREATE0(Class) \
  static typename Traits <Class>::Ptr_t create () { \
    typename Traits <Class>::Ptr_t ptr (new Class ()); \
    ptr->init (ptr); \
    return ptr; \
  }

static const int Nv = 9;
typedef Eigen::Matrix<value_type, Nv, 1> vector9_t;
typedef solver::lineSearch::Backtracking LineSearch_t;
using hpp::constraints::solver::HierarchicalIterative;
using boost::assign::list_of;
using boost::shared_ptr;
using hpp::setpyformat;
using hpp::one_line;

namespace hpp
{
namespace constraints
{

class PointA;
class PointB;
class PointC;
class PointD;

class PointB : public CalculusBase <PointB, vector3_t, JacobianMatrix>
{
  public:
    HPP_CONSTRAINTS_CB_CREATE0(PointB)
    DECLARE_EXPRESSION(PointB)

    PointB (const CalculusBase<PointB, vector3_t, JacobianMatrix>& other) :
      CalculusBase <PointB, eigen::vector3_t, JacobianMatrix> (other) {}

    PointB (const PointB& point) :
      CalculusBase <PointB, vector3_t, JacobianMatrix> (point) {}

    PointB () : CalculusBase <PointB, vector3_t, JacobianMatrix> ()
    {
      this->jacobian_.resize(3,Nv);
      this->jacobian_.rightCols<   3>().setIdentity();
      this->jacobian_.leftCols <Nv-3>().setZero();
    }

    void impl_value (const ConfigurationIn_t x)
    {
      this->value_ = x.head<3>();
    }
    void impl_jacobian (const ConfigurationIn_t ) {}
};

class PointC : public CalculusBase <PointC, vector3_t, JacobianMatrix>
{
  public:
    HPP_CONSTRAINTS_CB_CREATE0(PointC)
    DECLARE_EXPRESSION(PointC)

    PointC (const CalculusBase<PointC, vector3_t, JacobianMatrix>& other) :
      CalculusBase <PointC, eigen::vector3_t, JacobianMatrix> (other) {}

    PointC (const PointC& point) :
      CalculusBase <PointC, vector3_t, JacobianMatrix> (point) {}

    PointC () : CalculusBase <PointC, vector3_t, JacobianMatrix> ()
    {
      this->jacobian_.resize(3,Nv);
      this->jacobian_.rightCols<   3>().setIdentity();
      this->jacobian_.leftCols <Nv-3>().setZero();
      this->jacobian_(0, 3) = 1.;
    }

    void impl_value (const ConfigurationIn_t x)
    {
      this->value_ = x.head<3>();
      this->value_[0] += x[3];
    }
    void impl_jacobian (const ConfigurationIn_t ) {}
};

class PointD : public CalculusBase <PointD, vector3_t, JacobianMatrix>
{
  public:
    HPP_CONSTRAINTS_CB_CREATE0(PointD)
    DECLARE_EXPRESSION(PointD)

    PointD (const CalculusBase<PointD, vector3_t, JacobianMatrix>& other) :
      CalculusBase <PointD, eigen::vector3_t, JacobianMatrix> (other) {}

    PointD (const PointD& point) :
      CalculusBase <PointD, vector3_t, JacobianMatrix> (point) {}

    PointD () : CalculusBase <PointD, vector3_t, JacobianMatrix> ()
    {
      this->jacobian_.resize(3,Nv);
      this->jacobian_.setZero();
      this->jacobian_(0, 0) = 1.;
      this->jacobian_(0, 4) = 1.;
      this->jacobian_(1, 5) = 1.;
      this->jacobian_(2, 2) = 1.;
    }

    void impl_value (const ConfigurationIn_t x)
    {
      this->value_ << x[0] + x[4], x[5], x[2];
    }
    void impl_jacobian (const ConfigurationIn_t ) {}
};

class PointA : public CalculusBase <PointA, vector3_t, JacobianMatrix>
{
  public:
    HPP_CONSTRAINTS_CB_CREATE0(PointA)
    DECLARE_EXPRESSION(PointA)

    PointA (const CalculusBase<PointA, vector3_t, JacobianMatrix>& other) :
      CalculusBase <PointA, eigen::vector3_t, JacobianMatrix> (other) {}

    PointA (const PointA& point) :
      CalculusBase <PointA, vector3_t, JacobianMatrix> (point) {}

    PointA () : CalculusBase <PointA, vector3_t, JacobianMatrix> ()
    {
      this->jacobian_.resize(3,Nv);
      this->jacobian_.rightCols<Nv-3>().setIdentity();
      this->jacobian_.leftCols <   3>().setZero();
      this->jacobian_(2, 2) = 1.;
    }

    void impl_value (const ConfigurationIn_t x)
    {
      this->value_ = x.tail<3>();
      this->value_[2] += x[2];
    }
    void impl_jacobian (const ConfigurationIn_t ) {}
};

} // namespace constraints
} // namespace hpp

typedef Difference<PointB,PointA> AB_t;
typedef Difference<PointC,PointA> AC_t;
typedef Difference<PointD,PointA> AD_t;

typedef ScalarProduct<AB_t,PointA> AB_OA_t;
typedef ScalarProduct<AC_t,PointA> AC_OA_t;
typedef ScalarProduct<AD_t,PointA> AD_OA_t;

typedef CrossProduct<AB_t,AC_t> ABC_t;
typedef CrossProduct<AC_t,AD_t> ACD_t;
typedef CrossProduct<AD_t,AB_t> ADB_t;

typedef ScalarProduct<ABC_t,PointA> ABC_OA_t;
typedef ScalarProduct<ACD_t,PointA> ACD_OA_t;
typedef ScalarProduct<ADB_t,PointA> ADB_OA_t;

typedef CrossProduct<ABC_t,AB_t> ABC_AB_t;
typedef CrossProduct<ABC_t,AC_t> ABC_AC_t;
typedef CrossProduct<ACD_t,AC_t> ACD_AC_t;
typedef CrossProduct<ACD_t,AD_t> ACD_AD_t;
typedef CrossProduct<ADB_t,AD_t> ADB_AD_t;
typedef CrossProduct<ADB_t,AB_t> ADB_AB_t;

typedef ScalarProduct<ABC_AB_t,PointA> ABC_AB_OA_t;
typedef ScalarProduct<ABC_AC_t,PointA> ABC_AC_OA_t;
typedef ScalarProduct<ACD_AC_t,PointA> ACD_AC_OA_t;
typedef ScalarProduct<ACD_AD_t,PointA> ACD_AD_OA_t;
typedef ScalarProduct<ADB_AD_t,PointA> ADB_AD_OA_t;
typedef ScalarProduct<ADB_AB_t,PointA> ADB_AB_OA_t;

//ComparisonTypes_t negative (1, Inferior);

template <typename Expression>
ImplicitPtr_t implicitConstraint (const std::string& name,
    const shared_ptr<Expression>& expr, ComparisonTypes_t comp)
{
  typedef SymbolicFunction<Expression> SF_t;
  return Implicit::create (SF_t::create (name, Nv, Nv, 1, expr), comp);
}

struct Data
{
  PointA::Ptr_t A;
  PointB::Ptr_t B;
  PointC::Ptr_t C;
  PointD::Ptr_t D;

  shared_ptr<AB_t> AB;
  shared_ptr<AC_t> AC;
  shared_ptr<AD_t> AD;

  shared_ptr<AB_OA_t> AB_OA;
  shared_ptr<AC_OA_t> CB_OA;
  shared_ptr<AD_OA_t> DB_OA;

  shared_ptr<ABC_t> ABC;
  shared_ptr<ACD_t> ACD;
  shared_ptr<ADB_t> ADB;

  shared_ptr<ABC_OA_t> ABC_OA;
  shared_ptr<ACD_OA_t> ACD_OA;
  shared_ptr<ADB_OA_t> ADB_OA;

  shared_ptr<ABC_AB_t> ABC_AB;
  shared_ptr<ABC_AC_t> ABC_AC;
  shared_ptr<ACD_AC_t> ACD_AC;
  shared_ptr<ACD_AD_t> ACD_AD;
  shared_ptr<ADB_AD_t> ADB_AD;
  shared_ptr<ADB_AB_t> ADB_AB;

  shared_ptr<ABC_AB_OA_t> ABC_AB_OA;
  shared_ptr<ABC_AC_OA_t> ABC_AC_OA;
  shared_ptr<ACD_AC_OA_t> ACD_AC_OA;
  shared_ptr<ACD_AD_OA_t> ACD_AD_OA;
  shared_ptr<ADB_AD_OA_t> ADB_AD_OA;
  shared_ptr<ADB_AB_OA_t> ADB_AB_OA;

  ImplicitPtr_t c0_AB_OA, c0_CB_OA, c0_DB_OA,
                c1_AB_OA, c1_CB_OA, c1_DB_OA,
                c0_ABC_OA, c0_ACD_OA, c0_ADB_OA,
                c1_ABC_OA, c1_ACD_OA, c1_ADB_OA,
                c0_ABC_AB_OA, c0_ABC_AC_OA,
                c1_ABC_AB_OA, c1_ABC_AC_OA,
                c0_ACD_AC_OA, c0_ACD_AD_OA,
                c1_ACD_AC_OA, c1_ACD_AD_OA,
                c0_ADB_AD_OA, c0_ADB_AB_OA,
                c1_ADB_AD_OA, c1_ADB_AB_OA;
  ImplicitPtr_t     a[12];
  ImplicitPtr_t not_a[12];

  ImplicitPtr_t O_inside_1, O_inside_2;

  Data ()
  {
    A =(PointA::create());
    B =(PointB::create());
    C =(PointC::create());
    D =(PointD::create());

    AB =(B-A);
    AC =(C-A);
    AD =(D-A);
    AB_OA =(AB*A);
    CB_OA =(AC*A);
    DB_OA =(AD*A);

    ABC =(AB^AC);
    ACD =(AC^AD);
    ADB =(AD^AB);

    ABC_OA =(ABC*A);
    ACD_OA =(ACD*A);
    ADB_OA =(ADB*A);

    ABC_AB =(ABC^AB);
    ABC_AC =(ABC^AC);
    ACD_AC =(ACD^AC);
    ACD_AD =(ACD^AD);
    ADB_AD =(ADB^AD);
    ADB_AB =(ADB^AB);

    ABC_AB_OA =(ABC_AB*A);
    ABC_AC_OA =(ABC_AC*A);
    ACD_AC_OA =(ACD_AC*A);
    ACD_AD_OA =(ACD_AD*A);
    ADB_AD_OA =(ADB_AD*A);
    ADB_AB_OA =(ADB_AB*A);

    ComparisonTypes_t negative (1, Inferior);
    ComparisonTypes_t positive (1, Superior);

    c1_AB_OA     = implicitConstraint ("AB.AO >= 0"    , AB_OA    , negative);
    c1_CB_OA     = implicitConstraint ("CB.AO >= 0"    , CB_OA    , negative);
    c1_DB_OA     = implicitConstraint ("DB.AO >= 0"    , DB_OA    , negative);
    c1_ABC_OA    = implicitConstraint ("ABC.AO >= 0"   , ABC_OA   , negative);
    c1_ACD_OA    = implicitConstraint ("ACD.AO >= 0"   , ACD_OA   , negative);
    c1_ADB_OA    = implicitConstraint ("ADB.AO >= 0"   , ADB_OA   , negative);
    c1_ABC_AB_OA = implicitConstraint ("ABC^AB.AO >= 0", ABC_AB_OA, negative);
    c1_ABC_AC_OA = implicitConstraint ("ABC^AC.AO >= 0", ABC_AC_OA, negative);
    c1_ACD_AC_OA = implicitConstraint ("ACD^AC.AO >= 0", ACD_AC_OA, negative);
    c1_ACD_AD_OA = implicitConstraint ("ACD^AD.AO >= 0", ACD_AD_OA, negative);
    c1_ADB_AD_OA = implicitConstraint ("ADB^AD.AO >= 0", ADB_AD_OA, negative);
    c1_ADB_AB_OA = implicitConstraint ("ADB^AB.AO >= 0", ADB_AB_OA, negative);

    c0_AB_OA     = implicitConstraint ("AB.AO <= 0"    , AB_OA    , positive);
    c0_CB_OA     = implicitConstraint ("CB.AO <= 0"    , CB_OA    , positive);
    c0_DB_OA     = implicitConstraint ("DB.AO <= 0"    , DB_OA    , positive);
    c0_ABC_OA    = implicitConstraint ("ABC.AO <= 0"   , ABC_OA   , positive);
    c0_ACD_OA    = implicitConstraint ("ACD.AO <= 0"   , ACD_OA   , positive);
    c0_ADB_OA    = implicitConstraint ("ADB.AO <= 0"   , ADB_OA   , positive);
    c0_ABC_AB_OA = implicitConstraint ("ABC^AB.AO <= 0", ABC_AB_OA, positive);
    c0_ABC_AC_OA = implicitConstraint ("ABC^AC.AO <= 0", ABC_AC_OA, positive);
    c0_ACD_AC_OA = implicitConstraint ("ACD^AC.AO <= 0", ACD_AC_OA, positive);
    c0_ACD_AD_OA = implicitConstraint ("ACD^AD.AO <= 0", ACD_AD_OA, positive);
    c0_ADB_AD_OA = implicitConstraint ("ADB^AD.AO <= 0", ADB_AD_OA, positive);
    c0_ADB_AB_OA = implicitConstraint ("ADB^AB.AO <= 0", ADB_AB_OA, positive);

    shared_ptr<Point> z = Point::create(vector3_t(0,0,1),Nv);
    O_inside_1 = implicitConstraint ("BD^z.OB >= 0", ( (D-B)^z )*B, negative);
    O_inside_2 = implicitConstraint ("CD^z.OC <= 0", ( (D-C)^z )*C, positive);

    a[ 0] = c1_ABC_OA   ;
    a[ 1] = c1_ACD_OA   ;
    a[ 2] = c1_ADB_OA   ;
    a[ 3] = c1_ABC_AB_OA;
    a[ 4] = c1_ABC_AC_OA;
    a[ 5] = c1_ACD_AC_OA;
    a[ 6] = c1_ACD_AD_OA;
    a[ 7] = c1_ADB_AD_OA;
    a[ 8] = c1_ADB_AB_OA;
    a[ 9] = c1_AB_OA    ;
    a[10] = c1_CB_OA    ;
    a[11] = c1_DB_OA    ;

    not_a[ 0] = c0_ABC_OA   ;
    not_a[ 1] = c0_ACD_OA   ;
    not_a[ 2] = c0_ADB_OA   ;
    not_a[ 3] = c0_ABC_AB_OA;
    not_a[ 4] = c0_ABC_AC_OA;
    not_a[ 5] = c0_ACD_AC_OA;
    not_a[ 6] = c0_ACD_AD_OA;
    not_a[ 7] = c0_ADB_AD_OA;
    not_a[ 8] = c0_ADB_AB_OA;
    not_a[ 9] = c0_AB_OA    ;
    not_a[10] = c0_CB_OA    ;
    not_a[11] = c0_DB_OA    ;
  }
};

#define CHECK_MSG(cond,msg) if(!(cond)) std::cerr << __FILE__ ":" << __LINE__ << ": " << msg << ": failed" << std::endl;
#define CHECK_EQUAL(a,b) CHECK_MSG(a==b, #a " (" << a << ") == " #b " (" << b << ")")
#define CHECK(cond) CHECK_MSG(cond,#cond)

bool saturation(vectorIn_t q, vectorOut_t qSat, Eigen::VectorXi& saturation)
{
  static const value_type eps = 1e-1;
  bool ret = false;
  qSat = q;
  saturation.setZero();
  if (q[1] > 0) { // b1 <= 0
    saturation[1] = 1;
    qSat[1] = 0.;
    ret = true;
  }
  if (q[2] > -eps) { // z < 0
    saturation[2] = 1;
    qSat[2] = -eps;
    ret = true;
  }
  if (q[3] < eps) { // s0 > 0
    saturation[3] = -1;
    qSat[3] = eps;
    ret = true;
  }
  if (q[5] < eps) { // d1 > 0
    saturation[5] = -1;
    qSat[5] = eps;
    ret = true;
  }
  if (q[8] < eps) { // h > 0
    saturation[8] = -1;
    qSat[8] = eps;
    ret = true;
  }
  return ret;
}

value_type ineq_thr = 0.001;
size_type max_iter = 100;

HierarchicalIterative createSolver (Data& data)
{
  HierarchicalIterative solver (LiegroupSpace::Rn(Nv));
  solver.saturation (saturation);
  solver.inequalityThreshold (ineq_thr);
  solver.maxIterations (max_iter);
  solver.add (data.O_inside_1, 0);
  solver.add (data.O_inside_2, 0);
  return solver;
}

void check (Data& d)
{
  vector9_t x;
  x.setRandom();

  d.ABC->invalidate ();
  d.ABC->computeValue (x);
  CHECK_EQUAL(d.A->value()[0], x(6));
  CHECK_EQUAL(d.A->value()[1], x(7));
  CHECK_EQUAL(d.A->value()[2], x(8)+x[2]);
  CHECK_EQUAL(d.B->value()   , x.head<3>());

  HierarchicalIterative solver (createSolver(d));
  solver.add (d.c1_AB_OA, 1);

  x << 0, -1, -1, 1, 0, 1, 0, 0, 2;

  HierarchicalIterative::Status status = solver.solve<LineSearch_t> (x);
  CHECK_EQUAL(status,HierarchicalIterative::SUCCESS);
}

int checkCase (Data& d, int N, std::vector<int> conds)
{
  HierarchicalIterative solver (createSolver(d));
  for (std::size_t i = 0; i < conds.size(); ++i) {
    if (i>0) std::cout << '.';
    if (conds[i] > 0) {
      solver.add (d.    a[ conds[i]-1], 1);
    } else {
      std::cout << '!';
      solver.add (d.not_a[-conds[i]-1], 1);
    }
    std::cout << 'a' << abs(conds[i]);
  }

  vector9_t x, y;
  bool success = false;
  Eigen::VectorXi sat(Nv);
  for (int i = 0; i < N && !success; ++i) {
    y.setRandom();
    saturation (y, x, sat);
    HierarchicalIterative::Status status = solver.solve<LineSearch_t> (x);
    if (status == HierarchicalIterative::SUCCESS)
      success = true;
  }
  CHECK(!saturation(x, y, sat));

  std::cout << " : ";
  if (success) {
    d.A->invalidate(); d.A->computeValue(x);
    d.B->invalidate(); d.B->computeValue(x);
    d.C->invalidate(); d.C->computeValue(x);
    d.D->invalidate(); d.D->computeValue(x);
    std::cout << "found.\nX = " << setpyformat << one_line(x) << '\n';
    std::cout << '\n'
      << "A = np.array(" << one_line(d.A->value()) << ")\n"
      << "B = np.array(" << one_line(d.B->value()) << ")\n"
      << "C = np.array(" << one_line(d.C->value()) << ")\n"
      << "D = np.array(" << one_line(d.D->value()) << ")\n"
      ;
  } else {
    std::cout << "failed\n";
  }

  std::cout << '\n';
  const NumericalConstraints_t& ncs = solver.constraints();
  for (std::size_t i = 0; i < ncs.size(); ++i) {
    if (i>0) std::cout << " . ";
    std::cout << ncs[i]->function().name();
  }
  std::cout << '\n';

  return (success ? 0 : 1);
}

template <typename T>
T parse (char* arg)
{
  std::istringstream iss (arg);
  T val;
  if (iss >> val)
    return val;
  std::ostringstream oss;
  oss << "Could not parse " << arg << " as an integer.";
  throw std::invalid_argument (oss.str());
}

int main (int argc, char** argv)
{
  int max_num_shoot = 100;
  std::vector<int> case_, defaultCase (list_of(1)(4)(-5).convert_to_container<std::vector<int> >());

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--num-shoot") == 0)
      max_num_shoot = parse<int> (argv[++i]);
    else if (strcmp(argv[i], "--ineq-thr") == 0)
      ineq_thr = parse<value_type> (argv[++i]);
    else if (strcmp(argv[i], "--max-iter") == 0)
      max_iter = parse<size_type> (argv[++i]);
    else
      case_.push_back(parse<int>(argv[i]));
  }
  if (case_.empty()) case_ = defaultCase;

  Data data;
  check (data);

  return checkCase (data, max_num_shoot, case_);
}
