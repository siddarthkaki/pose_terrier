// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: pose.proto

#include "pose.pb.h"

#include <algorithm>

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/stubs/port.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/wire_format_lite_inl.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// This is a temporary google only hack
#ifdef GOOGLE_PROTOBUF_ENFORCE_UNIQUENESS
#include "third_party/protobuf/version.h"
#endif
// @@protoc_insertion_point(includes)

namespace protobuf_pose_2eproto {
extern PROTOBUF_INTERNAL_EXPORT_protobuf_pose_2eproto ::google::protobuf::internal::SCCInfo<0> scc_info_Attitude;
extern PROTOBUF_INTERNAL_EXPORT_protobuf_pose_2eproto ::google::protobuf::internal::SCCInfo<0> scc_info_Position;
}  // namespace protobuf_pose_2eproto
namespace ProtoPose {
class PositionDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<Position>
      _instance;
} _Position_default_instance_;
class AttitudeDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<Attitude>
      _instance;
} _Attitude_default_instance_;
class PoseDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<Pose>
      _instance;
} _Pose_default_instance_;
}  // namespace ProtoPose
namespace protobuf_pose_2eproto {
static void InitDefaultsPosition() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ProtoPose::_Position_default_instance_;
    new (ptr) ::ProtoPose::Position();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ProtoPose::Position::InitAsDefaultInstance();
}

::google::protobuf::internal::SCCInfo<0> scc_info_Position =
    {{ATOMIC_VAR_INIT(::google::protobuf::internal::SCCInfoBase::kUninitialized), 0, InitDefaultsPosition}, {}};

static void InitDefaultsAttitude() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ProtoPose::_Attitude_default_instance_;
    new (ptr) ::ProtoPose::Attitude();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ProtoPose::Attitude::InitAsDefaultInstance();
}

::google::protobuf::internal::SCCInfo<0> scc_info_Attitude =
    {{ATOMIC_VAR_INIT(::google::protobuf::internal::SCCInfoBase::kUninitialized), 0, InitDefaultsAttitude}, {}};

static void InitDefaultsPose() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ProtoPose::_Pose_default_instance_;
    new (ptr) ::ProtoPose::Pose();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ProtoPose::Pose::InitAsDefaultInstance();
}

::google::protobuf::internal::SCCInfo<2> scc_info_Pose =
    {{ATOMIC_VAR_INIT(::google::protobuf::internal::SCCInfoBase::kUninitialized), 2, InitDefaultsPose}, {
      &protobuf_pose_2eproto::scc_info_Position.base,
      &protobuf_pose_2eproto::scc_info_Attitude.base,}};

void InitDefaults() {
  ::google::protobuf::internal::InitSCC(&scc_info_Position.base);
  ::google::protobuf::internal::InitSCC(&scc_info_Attitude.base);
  ::google::protobuf::internal::InitSCC(&scc_info_Pose.base);
}

::google::protobuf::Metadata file_level_metadata[3];

const ::google::protobuf::uint32 TableStruct::offsets[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Position, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Position, x_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Position, y_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Position, z_),
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Attitude, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Attitude, qw_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Attitude, qx_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Attitude, qy_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Attitude, qz_),
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Pose, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Pose, pos_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Pose, att_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoPose::Pose, time_stamp_),
};
static const ::google::protobuf::internal::MigrationSchema schemas[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, sizeof(::ProtoPose::Position)},
  { 8, -1, sizeof(::ProtoPose::Attitude)},
  { 17, -1, sizeof(::ProtoPose::Pose)},
};

static ::google::protobuf::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::google::protobuf::Message*>(&::ProtoPose::_Position_default_instance_),
  reinterpret_cast<const ::google::protobuf::Message*>(&::ProtoPose::_Attitude_default_instance_),
  reinterpret_cast<const ::google::protobuf::Message*>(&::ProtoPose::_Pose_default_instance_),
};

void protobuf_AssignDescriptors() {
  AddDescriptors();
  AssignDescriptors(
      "pose.proto", schemas, file_default_instances, TableStruct::offsets,
      file_level_metadata, NULL, NULL);
}

void protobuf_AssignDescriptorsOnce() {
  static ::google::protobuf::internal::once_flag once;
  ::google::protobuf::internal::call_once(once, protobuf_AssignDescriptors);
}

void protobuf_RegisterTypes(const ::std::string&) GOOGLE_PROTOBUF_ATTRIBUTE_COLD;
void protobuf_RegisterTypes(const ::std::string&) {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::internal::RegisterAllTypes(file_level_metadata, 3);
}

void AddDescriptorsImpl() {
  InitDefaults();
  static const char descriptor[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
      "\n\npose.proto\022\tProtoPose\"+\n\010Position\022\t\n\001x"
      "\030\001 \001(\001\022\t\n\001y\030\002 \001(\001\022\t\n\001z\030\003 \001(\001\":\n\010Attitude"
      "\022\n\n\002qw\030\001 \001(\001\022\n\n\002qx\030\002 \001(\001\022\n\n\002qy\030\003 \001(\001\022\n\n\002"
      "qz\030\004 \001(\001\"^\n\004Pose\022 \n\003pos\030\001 \001(\0132\023.ProtoPos"
      "e.Position\022 \n\003att\030\002 \001(\0132\023.ProtoPose.Atti"
      "tude\022\022\n\ntime_stamp\030\003 \001(\001b\006proto3"
  };
  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
      descriptor, 232);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "pose.proto", &protobuf_RegisterTypes);
}

void AddDescriptors() {
  static ::google::protobuf::internal::once_flag once;
  ::google::protobuf::internal::call_once(once, AddDescriptorsImpl);
}
// Force AddDescriptors() to be called at dynamic initialization time.
struct StaticDescriptorInitializer {
  StaticDescriptorInitializer() {
    AddDescriptors();
  }
} static_descriptor_initializer;
}  // namespace protobuf_pose_2eproto
namespace ProtoPose {

// ===================================================================

void Position::InitAsDefaultInstance() {
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Position::kXFieldNumber;
const int Position::kYFieldNumber;
const int Position::kZFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Position::Position()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  ::google::protobuf::internal::InitSCC(
      &protobuf_pose_2eproto::scc_info_Position.base);
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoPose.Position)
}
Position::Position(const Position& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::memcpy(&x_, &from.x_,
    static_cast<size_t>(reinterpret_cast<char*>(&z_) -
    reinterpret_cast<char*>(&x_)) + sizeof(z_));
  // @@protoc_insertion_point(copy_constructor:ProtoPose.Position)
}

void Position::SharedCtor() {
  ::memset(&x_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&z_) -
      reinterpret_cast<char*>(&x_)) + sizeof(z_));
}

Position::~Position() {
  // @@protoc_insertion_point(destructor:ProtoPose.Position)
  SharedDtor();
}

void Position::SharedDtor() {
}

void Position::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ::google::protobuf::Descriptor* Position::descriptor() {
  ::protobuf_pose_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_pose_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const Position& Position::default_instance() {
  ::google::protobuf::internal::InitSCC(&protobuf_pose_2eproto::scc_info_Position.base);
  return *internal_default_instance();
}


void Position::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoPose.Position)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  ::memset(&x_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&z_) -
      reinterpret_cast<char*>(&x_)) + sizeof(z_));
  _internal_metadata_.Clear();
}

bool Position::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoPose.Position)
  for (;;) {
    ::std::pair<::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // double x = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(9u /* 9 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &x_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double y = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(17u /* 17 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &y_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double z = 3;
      case 3: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(25u /* 25 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &z_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, _internal_metadata_.mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:ProtoPose.Position)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:ProtoPose.Position)
  return false;
#undef DO_
}

void Position::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:ProtoPose.Position)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double x = 1;
  if (this->x() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(1, this->x(), output);
  }

  // double y = 2;
  if (this->y() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(2, this->y(), output);
  }

  // double z = 3;
  if (this->z() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(3, this->z(), output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:ProtoPose.Position)
}

::google::protobuf::uint8* Position::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:ProtoPose.Position)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double x = 1;
  if (this->x() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(1, this->x(), target);
  }

  // double y = 2;
  if (this->y() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(2, this->y(), target);
  }

  // double z = 3;
  if (this->z() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(3, this->z(), target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ProtoPose.Position)
  return target;
}

size_t Position::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoPose.Position)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // double x = 1;
  if (this->x() != 0) {
    total_size += 1 + 8;
  }

  // double y = 2;
  if (this->y() != 0) {
    total_size += 1 + 8;
  }

  // double z = 3;
  if (this->z() != 0) {
    total_size += 1 + 8;
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void Position::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoPose.Position)
  GOOGLE_DCHECK_NE(&from, this);
  const Position* source =
      ::google::protobuf::internal::DynamicCastToGenerated<const Position>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ProtoPose.Position)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ProtoPose.Position)
    MergeFrom(*source);
  }
}

void Position::MergeFrom(const Position& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ProtoPose.Position)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  if (from.x() != 0) {
    set_x(from.x());
  }
  if (from.y() != 0) {
    set_y(from.y());
  }
  if (from.z() != 0) {
    set_z(from.z());
  }
}

void Position::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ProtoPose.Position)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Position::CopyFrom(const Position& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ProtoPose.Position)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Position::IsInitialized() const {
  return true;
}

void Position::Swap(Position* other) {
  if (other == this) return;
  InternalSwap(other);
}
void Position::InternalSwap(Position* other) {
  using std::swap;
  swap(x_, other->x_);
  swap(y_, other->y_);
  swap(z_, other->z_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
}

::google::protobuf::Metadata Position::GetMetadata() const {
  protobuf_pose_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_pose_2eproto::file_level_metadata[kIndexInFileMessages];
}


// ===================================================================

void Attitude::InitAsDefaultInstance() {
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Attitude::kQwFieldNumber;
const int Attitude::kQxFieldNumber;
const int Attitude::kQyFieldNumber;
const int Attitude::kQzFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Attitude::Attitude()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  ::google::protobuf::internal::InitSCC(
      &protobuf_pose_2eproto::scc_info_Attitude.base);
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoPose.Attitude)
}
Attitude::Attitude(const Attitude& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::memcpy(&qw_, &from.qw_,
    static_cast<size_t>(reinterpret_cast<char*>(&qz_) -
    reinterpret_cast<char*>(&qw_)) + sizeof(qz_));
  // @@protoc_insertion_point(copy_constructor:ProtoPose.Attitude)
}

void Attitude::SharedCtor() {
  ::memset(&qw_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&qz_) -
      reinterpret_cast<char*>(&qw_)) + sizeof(qz_));
}

Attitude::~Attitude() {
  // @@protoc_insertion_point(destructor:ProtoPose.Attitude)
  SharedDtor();
}

void Attitude::SharedDtor() {
}

void Attitude::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ::google::protobuf::Descriptor* Attitude::descriptor() {
  ::protobuf_pose_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_pose_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const Attitude& Attitude::default_instance() {
  ::google::protobuf::internal::InitSCC(&protobuf_pose_2eproto::scc_info_Attitude.base);
  return *internal_default_instance();
}


void Attitude::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoPose.Attitude)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  ::memset(&qw_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&qz_) -
      reinterpret_cast<char*>(&qw_)) + sizeof(qz_));
  _internal_metadata_.Clear();
}

bool Attitude::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoPose.Attitude)
  for (;;) {
    ::std::pair<::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // double qw = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(9u /* 9 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &qw_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double qx = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(17u /* 17 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &qx_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double qy = 3;
      case 3: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(25u /* 25 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &qy_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double qz = 4;
      case 4: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(33u /* 33 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &qz_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, _internal_metadata_.mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:ProtoPose.Attitude)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:ProtoPose.Attitude)
  return false;
#undef DO_
}

void Attitude::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:ProtoPose.Attitude)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double qw = 1;
  if (this->qw() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(1, this->qw(), output);
  }

  // double qx = 2;
  if (this->qx() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(2, this->qx(), output);
  }

  // double qy = 3;
  if (this->qy() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(3, this->qy(), output);
  }

  // double qz = 4;
  if (this->qz() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(4, this->qz(), output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:ProtoPose.Attitude)
}

::google::protobuf::uint8* Attitude::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:ProtoPose.Attitude)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double qw = 1;
  if (this->qw() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(1, this->qw(), target);
  }

  // double qx = 2;
  if (this->qx() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(2, this->qx(), target);
  }

  // double qy = 3;
  if (this->qy() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(3, this->qy(), target);
  }

  // double qz = 4;
  if (this->qz() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(4, this->qz(), target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ProtoPose.Attitude)
  return target;
}

size_t Attitude::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoPose.Attitude)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // double qw = 1;
  if (this->qw() != 0) {
    total_size += 1 + 8;
  }

  // double qx = 2;
  if (this->qx() != 0) {
    total_size += 1 + 8;
  }

  // double qy = 3;
  if (this->qy() != 0) {
    total_size += 1 + 8;
  }

  // double qz = 4;
  if (this->qz() != 0) {
    total_size += 1 + 8;
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void Attitude::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoPose.Attitude)
  GOOGLE_DCHECK_NE(&from, this);
  const Attitude* source =
      ::google::protobuf::internal::DynamicCastToGenerated<const Attitude>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ProtoPose.Attitude)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ProtoPose.Attitude)
    MergeFrom(*source);
  }
}

void Attitude::MergeFrom(const Attitude& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ProtoPose.Attitude)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  if (from.qw() != 0) {
    set_qw(from.qw());
  }
  if (from.qx() != 0) {
    set_qx(from.qx());
  }
  if (from.qy() != 0) {
    set_qy(from.qy());
  }
  if (from.qz() != 0) {
    set_qz(from.qz());
  }
}

void Attitude::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ProtoPose.Attitude)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Attitude::CopyFrom(const Attitude& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ProtoPose.Attitude)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Attitude::IsInitialized() const {
  return true;
}

void Attitude::Swap(Attitude* other) {
  if (other == this) return;
  InternalSwap(other);
}
void Attitude::InternalSwap(Attitude* other) {
  using std::swap;
  swap(qw_, other->qw_);
  swap(qx_, other->qx_);
  swap(qy_, other->qy_);
  swap(qz_, other->qz_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
}

::google::protobuf::Metadata Attitude::GetMetadata() const {
  protobuf_pose_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_pose_2eproto::file_level_metadata[kIndexInFileMessages];
}


// ===================================================================

void Pose::InitAsDefaultInstance() {
  ::ProtoPose::_Pose_default_instance_._instance.get_mutable()->pos_ = const_cast< ::ProtoPose::Position*>(
      ::ProtoPose::Position::internal_default_instance());
  ::ProtoPose::_Pose_default_instance_._instance.get_mutable()->att_ = const_cast< ::ProtoPose::Attitude*>(
      ::ProtoPose::Attitude::internal_default_instance());
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Pose::kPosFieldNumber;
const int Pose::kAttFieldNumber;
const int Pose::kTimeStampFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Pose::Pose()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  ::google::protobuf::internal::InitSCC(
      &protobuf_pose_2eproto::scc_info_Pose.base);
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoPose.Pose)
}
Pose::Pose(const Pose& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  if (from.has_pos()) {
    pos_ = new ::ProtoPose::Position(*from.pos_);
  } else {
    pos_ = NULL;
  }
  if (from.has_att()) {
    att_ = new ::ProtoPose::Attitude(*from.att_);
  } else {
    att_ = NULL;
  }
  time_stamp_ = from.time_stamp_;
  // @@protoc_insertion_point(copy_constructor:ProtoPose.Pose)
}

void Pose::SharedCtor() {
  ::memset(&pos_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&time_stamp_) -
      reinterpret_cast<char*>(&pos_)) + sizeof(time_stamp_));
}

Pose::~Pose() {
  // @@protoc_insertion_point(destructor:ProtoPose.Pose)
  SharedDtor();
}

void Pose::SharedDtor() {
  if (this != internal_default_instance()) delete pos_;
  if (this != internal_default_instance()) delete att_;
}

void Pose::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ::google::protobuf::Descriptor* Pose::descriptor() {
  ::protobuf_pose_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_pose_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const Pose& Pose::default_instance() {
  ::google::protobuf::internal::InitSCC(&protobuf_pose_2eproto::scc_info_Pose.base);
  return *internal_default_instance();
}


void Pose::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoPose.Pose)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  if (GetArenaNoVirtual() == NULL && pos_ != NULL) {
    delete pos_;
  }
  pos_ = NULL;
  if (GetArenaNoVirtual() == NULL && att_ != NULL) {
    delete att_;
  }
  att_ = NULL;
  time_stamp_ = 0;
  _internal_metadata_.Clear();
}

bool Pose::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoPose.Pose)
  for (;;) {
    ::std::pair<::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // .ProtoPose.Position pos = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(10u /* 10 & 0xFF */)) {
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessage(
               input, mutable_pos()));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // .ProtoPose.Attitude att = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(18u /* 18 & 0xFF */)) {
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessage(
               input, mutable_att()));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double time_stamp = 3;
      case 3: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(25u /* 25 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &time_stamp_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, _internal_metadata_.mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:ProtoPose.Pose)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:ProtoPose.Pose)
  return false;
#undef DO_
}

void Pose::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:ProtoPose.Pose)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // .ProtoPose.Position pos = 1;
  if (this->has_pos()) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      1, this->_internal_pos(), output);
  }

  // .ProtoPose.Attitude att = 2;
  if (this->has_att()) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      2, this->_internal_att(), output);
  }

  // double time_stamp = 3;
  if (this->time_stamp() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(3, this->time_stamp(), output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:ProtoPose.Pose)
}

::google::protobuf::uint8* Pose::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:ProtoPose.Pose)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // .ProtoPose.Position pos = 1;
  if (this->has_pos()) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageToArray(
        1, this->_internal_pos(), deterministic, target);
  }

  // .ProtoPose.Attitude att = 2;
  if (this->has_att()) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageToArray(
        2, this->_internal_att(), deterministic, target);
  }

  // double time_stamp = 3;
  if (this->time_stamp() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(3, this->time_stamp(), target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ProtoPose.Pose)
  return target;
}

size_t Pose::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoPose.Pose)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // .ProtoPose.Position pos = 1;
  if (this->has_pos()) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::MessageSize(
        *pos_);
  }

  // .ProtoPose.Attitude att = 2;
  if (this->has_att()) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::MessageSize(
        *att_);
  }

  // double time_stamp = 3;
  if (this->time_stamp() != 0) {
    total_size += 1 + 8;
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void Pose::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoPose.Pose)
  GOOGLE_DCHECK_NE(&from, this);
  const Pose* source =
      ::google::protobuf::internal::DynamicCastToGenerated<const Pose>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ProtoPose.Pose)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ProtoPose.Pose)
    MergeFrom(*source);
  }
}

void Pose::MergeFrom(const Pose& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ProtoPose.Pose)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  if (from.has_pos()) {
    mutable_pos()->::ProtoPose::Position::MergeFrom(from.pos());
  }
  if (from.has_att()) {
    mutable_att()->::ProtoPose::Attitude::MergeFrom(from.att());
  }
  if (from.time_stamp() != 0) {
    set_time_stamp(from.time_stamp());
  }
}

void Pose::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ProtoPose.Pose)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Pose::CopyFrom(const Pose& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ProtoPose.Pose)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Pose::IsInitialized() const {
  return true;
}

void Pose::Swap(Pose* other) {
  if (other == this) return;
  InternalSwap(other);
}
void Pose::InternalSwap(Pose* other) {
  using std::swap;
  swap(pos_, other->pos_);
  swap(att_, other->att_);
  swap(time_stamp_, other->time_stamp_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
}

::google::protobuf::Metadata Pose::GetMetadata() const {
  protobuf_pose_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_pose_2eproto::file_level_metadata[kIndexInFileMessages];
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace ProtoPose
namespace google {
namespace protobuf {
template<> GOOGLE_PROTOBUF_ATTRIBUTE_NOINLINE ::ProtoPose::Position* Arena::CreateMaybeMessage< ::ProtoPose::Position >(Arena* arena) {
  return Arena::CreateInternal< ::ProtoPose::Position >(arena);
}
template<> GOOGLE_PROTOBUF_ATTRIBUTE_NOINLINE ::ProtoPose::Attitude* Arena::CreateMaybeMessage< ::ProtoPose::Attitude >(Arena* arena) {
  return Arena::CreateInternal< ::ProtoPose::Attitude >(arena);
}
template<> GOOGLE_PROTOBUF_ATTRIBUTE_NOINLINE ::ProtoPose::Pose* Arena::CreateMaybeMessage< ::ProtoPose::Pose >(Arena* arena) {
  return Arena::CreateInternal< ::ProtoPose::Pose >(arena);
}
}  // namespace protobuf
}  // namespace google

// @@protoc_insertion_point(global_scope)
