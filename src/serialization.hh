#ifndef SRC_SERIALIZATION
#define SRC_SERIALIZATION

namespace hpp {
namespace constraints {
namespace internal {
template <typename Archive>
inline void serialize_joint(Archive& ar, const char* name, JointConstPtr_t& j) {
  JointPtr_t joint;
  if (Archive::is_saving::value)
    joint = const_pointer_cast<hpp::pinocchio::Joint>(j);
  ar& boost::serialization::make_nvp(name, joint);
  if (!Archive::is_saving::value) j = joint;
}
}  // namespace internal
}  // namespace constraints
}  // namespace hpp

#endif  // SRC_SERIALIZATION
