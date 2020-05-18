#ifndef SRC_SERIALIZATION
#define SRC_SERIALIZATION

namespace hpp {
namespace constraints {
namespace internal {
template<typename Archive>
inline void serialize_joint (Archive & ar, const char* name, JointConstPtr_t j)
{
  JointPtr_t joint = boost::const_pointer_cast<hpp::pinocchio::Joint>(j);
  ar & boost::serialization::make_nvp(name, joint);
  j = joint;
}
}
}
}

#endif // SRC_SERIALIZATION
