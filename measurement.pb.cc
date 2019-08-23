// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: measurement.proto

#define INTERNAL_SUPPRESS_PROTOBUF_FIELD_DEPRECATION
#include "measurement.pb.h"

#include <algorithm>

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/stubs/port.h>
#include <google/protobuf/stubs/once.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/wire_format_lite_inl.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)

namespace ProtoMeas {

namespace {

const ::google::protobuf::Descriptor* Position_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  Position_reflection_ = NULL;
const ::google::protobuf::Descriptor* Bearing_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  Bearing_reflection_ = NULL;
const ::google::protobuf::Descriptor* Measurements_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  Measurements_reflection_ = NULL;

}  // namespace


void protobuf_AssignDesc_measurement_2eproto() GOOGLE_ATTRIBUTE_COLD;
void protobuf_AssignDesc_measurement_2eproto() {
  protobuf_AddDesc_measurement_2eproto();
  const ::google::protobuf::FileDescriptor* file =
    ::google::protobuf::DescriptorPool::generated_pool()->FindFileByName(
      "measurement.proto");
  GOOGLE_CHECK(file != NULL);
  Position_descriptor_ = file->message_type(0);
  static const int Position_offsets_[3] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Position, x_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Position, y_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Position, z_),
  };
  Position_reflection_ =
    ::google::protobuf::internal::GeneratedMessageReflection::NewGeneratedMessageReflection(
      Position_descriptor_,
      Position::default_instance_,
      Position_offsets_,
      -1,
      -1,
      -1,
      sizeof(Position),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Position, _internal_metadata_),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Position, _is_default_instance_));
  Bearing_descriptor_ = file->message_type(1);
  static const int Bearing_offsets_[2] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Bearing, az_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Bearing, el_),
  };
  Bearing_reflection_ =
    ::google::protobuf::internal::GeneratedMessageReflection::NewGeneratedMessageReflection(
      Bearing_descriptor_,
      Bearing::default_instance_,
      Bearing_offsets_,
      -1,
      -1,
      -1,
      sizeof(Bearing),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Bearing, _internal_metadata_),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Bearing, _is_default_instance_));
  Measurements_descriptor_ = file->message_type(2);
  static const int Measurements_offsets_[3] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Measurements, num_feature_points_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Measurements, feature_points_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Measurements, bearings_),
  };
  Measurements_reflection_ =
    ::google::protobuf::internal::GeneratedMessageReflection::NewGeneratedMessageReflection(
      Measurements_descriptor_,
      Measurements::default_instance_,
      Measurements_offsets_,
      -1,
      -1,
      -1,
      sizeof(Measurements),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Measurements, _internal_metadata_),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Measurements, _is_default_instance_));
}

namespace {

GOOGLE_PROTOBUF_DECLARE_ONCE(protobuf_AssignDescriptors_once_);
inline void protobuf_AssignDescriptorsOnce() {
  ::google::protobuf::GoogleOnceInit(&protobuf_AssignDescriptors_once_,
                 &protobuf_AssignDesc_measurement_2eproto);
}

void protobuf_RegisterTypes(const ::std::string&) GOOGLE_ATTRIBUTE_COLD;
void protobuf_RegisterTypes(const ::std::string&) {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
      Position_descriptor_, &Position::default_instance());
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
      Bearing_descriptor_, &Bearing::default_instance());
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
      Measurements_descriptor_, &Measurements::default_instance());
}

}  // namespace

void protobuf_ShutdownFile_measurement_2eproto() {
  delete Position::default_instance_;
  delete Position_reflection_;
  delete Bearing::default_instance_;
  delete Bearing_reflection_;
  delete Measurements::default_instance_;
  delete Measurements_reflection_;
}

void protobuf_AddDesc_measurement_2eproto() GOOGLE_ATTRIBUTE_COLD;
void protobuf_AddDesc_measurement_2eproto() {
  static bool already_here = false;
  if (already_here) return;
  already_here = true;
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
    "\n\021measurement.proto\022\tProtoMeas\"+\n\010Positi"
    "on\022\t\n\001x\030\001 \001(\001\022\t\n\001y\030\002 \001(\001\022\t\n\001z\030\003 \001(\001\"!\n\007B"
    "earing\022\n\n\002az\030\001 \001(\001\022\n\n\002el\030\002 \001(\001\"}\n\014Measur"
    "ements\022\032\n\022num_feature_points\030\001 \001(\r\022+\n\016fe"
    "ature_points\030\002 \003(\0132\023.ProtoMeas.Position\022"
    "$\n\010bearings\030\003 \003(\0132\022.ProtoMeas.Bearingb\006p"
    "roto3", 245);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "measurement.proto", &protobuf_RegisterTypes);
  Position::default_instance_ = new Position();
  Bearing::default_instance_ = new Bearing();
  Measurements::default_instance_ = new Measurements();
  Position::default_instance_->InitAsDefaultInstance();
  Bearing::default_instance_->InitAsDefaultInstance();
  Measurements::default_instance_->InitAsDefaultInstance();
  ::google::protobuf::internal::OnShutdown(&protobuf_ShutdownFile_measurement_2eproto);
}

// Force AddDescriptors() to be called at static initialization time.
struct StaticDescriptorInitializer_measurement_2eproto {
  StaticDescriptorInitializer_measurement_2eproto() {
    protobuf_AddDesc_measurement_2eproto();
  }
} static_descriptor_initializer_measurement_2eproto_;

// ===================================================================

#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Position::kXFieldNumber;
const int Position::kYFieldNumber;
const int Position::kZFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Position::Position()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoMeas.Position)
}

void Position::InitAsDefaultInstance() {
  _is_default_instance_ = true;
}

Position::Position(const Position& from)
  : ::google::protobuf::Message(),
    _internal_metadata_(NULL) {
  SharedCtor();
  MergeFrom(from);
  // @@protoc_insertion_point(copy_constructor:ProtoMeas.Position)
}

void Position::SharedCtor() {
    _is_default_instance_ = false;
  _cached_size_ = 0;
  x_ = 0;
  y_ = 0;
  z_ = 0;
}

Position::~Position() {
  // @@protoc_insertion_point(destructor:ProtoMeas.Position)
  SharedDtor();
}

void Position::SharedDtor() {
  if (this != default_instance_) {
  }
}

void Position::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* Position::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return Position_descriptor_;
}

const Position& Position::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_measurement_2eproto();
  return *default_instance_;
}

Position* Position::default_instance_ = NULL;

Position* Position::New(::google::protobuf::Arena* arena) const {
  Position* n = new Position;
  if (arena != NULL) {
    arena->Own(n);
  }
  return n;
}

void Position::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoMeas.Position)
#if defined(__clang__)
#define ZR_HELPER_(f) \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Winvalid-offsetof\"") \
  __builtin_offsetof(Position, f) \
  _Pragma("clang diagnostic pop")
#else
#define ZR_HELPER_(f) reinterpret_cast<char*>(\
  &reinterpret_cast<Position*>(16)->f)
#endif

#define ZR_(first, last) do {\
  ::memset(&first, 0,\
           ZR_HELPER_(last) - ZR_HELPER_(first) + sizeof(last));\
} while (0)

  ZR_(x_, z_);

#undef ZR_HELPER_
#undef ZR_

}

bool Position::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoMeas.Position)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoff(127);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // optional double x = 1;
      case 1: {
        if (tag == 9) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &x_)));

        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(17)) goto parse_y;
        break;
      }

      // optional double y = 2;
      case 2: {
        if (tag == 17) {
         parse_y:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &y_)));

        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(25)) goto parse_z;
        break;
      }

      // optional double z = 3;
      case 3: {
        if (tag == 25) {
         parse_z:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &z_)));

        } else {
          goto handle_unusual;
        }
        if (input->ExpectAtEnd()) goto success;
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0 ||
            ::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormatLite::SkipField(input, tag));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:ProtoMeas.Position)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:ProtoMeas.Position)
  return false;
#undef DO_
}

void Position::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:ProtoMeas.Position)
  // optional double x = 1;
  if (this->x() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(1, this->x(), output);
  }

  // optional double y = 2;
  if (this->y() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(2, this->y(), output);
  }

  // optional double z = 3;
  if (this->z() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(3, this->z(), output);
  }

  // @@protoc_insertion_point(serialize_end:ProtoMeas.Position)
}

::google::protobuf::uint8* Position::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  // @@protoc_insertion_point(serialize_to_array_start:ProtoMeas.Position)
  // optional double x = 1;
  if (this->x() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(1, this->x(), target);
  }

  // optional double y = 2;
  if (this->y() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(2, this->y(), target);
  }

  // optional double z = 3;
  if (this->z() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(3, this->z(), target);
  }

  // @@protoc_insertion_point(serialize_to_array_end:ProtoMeas.Position)
  return target;
}

int Position::ByteSize() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoMeas.Position)
  int total_size = 0;

  // optional double x = 1;
  if (this->x() != 0) {
    total_size += 1 + 8;
  }

  // optional double y = 2;
  if (this->y() != 0) {
    total_size += 1 + 8;
  }

  // optional double z = 3;
  if (this->z() != 0) {
    total_size += 1 + 8;
  }

  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void Position::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoMeas.Position)
  if (GOOGLE_PREDICT_FALSE(&from == this)) {
    ::google::protobuf::internal::MergeFromFail(__FILE__, __LINE__);
  }
  const Position* source = 
      ::google::protobuf::internal::DynamicCastToGenerated<const Position>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ProtoMeas.Position)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ProtoMeas.Position)
    MergeFrom(*source);
  }
}

void Position::MergeFrom(const Position& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ProtoMeas.Position)
  if (GOOGLE_PREDICT_FALSE(&from == this)) {
    ::google::protobuf::internal::MergeFromFail(__FILE__, __LINE__);
  }
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
// @@protoc_insertion_point(generalized_copy_from_start:ProtoMeas.Position)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Position::CopyFrom(const Position& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ProtoMeas.Position)
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
  std::swap(x_, other->x_);
  std::swap(y_, other->y_);
  std::swap(z_, other->z_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
  std::swap(_cached_size_, other->_cached_size_);
}

::google::protobuf::Metadata Position::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = Position_descriptor_;
  metadata.reflection = Position_reflection_;
  return metadata;
}

#if PROTOBUF_INLINE_NOT_IN_HEADERS
// Position

// optional double x = 1;
void Position::clear_x() {
  x_ = 0;
}
 double Position::x() const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Position.x)
  return x_;
}
 void Position::set_x(double value) {
  
  x_ = value;
  // @@protoc_insertion_point(field_set:ProtoMeas.Position.x)
}

// optional double y = 2;
void Position::clear_y() {
  y_ = 0;
}
 double Position::y() const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Position.y)
  return y_;
}
 void Position::set_y(double value) {
  
  y_ = value;
  // @@protoc_insertion_point(field_set:ProtoMeas.Position.y)
}

// optional double z = 3;
void Position::clear_z() {
  z_ = 0;
}
 double Position::z() const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Position.z)
  return z_;
}
 void Position::set_z(double value) {
  
  z_ = value;
  // @@protoc_insertion_point(field_set:ProtoMeas.Position.z)
}

#endif  // PROTOBUF_INLINE_NOT_IN_HEADERS

// ===================================================================

#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Bearing::kAzFieldNumber;
const int Bearing::kElFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Bearing::Bearing()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoMeas.Bearing)
}

void Bearing::InitAsDefaultInstance() {
  _is_default_instance_ = true;
}

Bearing::Bearing(const Bearing& from)
  : ::google::protobuf::Message(),
    _internal_metadata_(NULL) {
  SharedCtor();
  MergeFrom(from);
  // @@protoc_insertion_point(copy_constructor:ProtoMeas.Bearing)
}

void Bearing::SharedCtor() {
    _is_default_instance_ = false;
  _cached_size_ = 0;
  az_ = 0;
  el_ = 0;
}

Bearing::~Bearing() {
  // @@protoc_insertion_point(destructor:ProtoMeas.Bearing)
  SharedDtor();
}

void Bearing::SharedDtor() {
  if (this != default_instance_) {
  }
}

void Bearing::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* Bearing::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return Bearing_descriptor_;
}

const Bearing& Bearing::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_measurement_2eproto();
  return *default_instance_;
}

Bearing* Bearing::default_instance_ = NULL;

Bearing* Bearing::New(::google::protobuf::Arena* arena) const {
  Bearing* n = new Bearing;
  if (arena != NULL) {
    arena->Own(n);
  }
  return n;
}

void Bearing::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoMeas.Bearing)
#if defined(__clang__)
#define ZR_HELPER_(f) \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Winvalid-offsetof\"") \
  __builtin_offsetof(Bearing, f) \
  _Pragma("clang diagnostic pop")
#else
#define ZR_HELPER_(f) reinterpret_cast<char*>(\
  &reinterpret_cast<Bearing*>(16)->f)
#endif

#define ZR_(first, last) do {\
  ::memset(&first, 0,\
           ZR_HELPER_(last) - ZR_HELPER_(first) + sizeof(last));\
} while (0)

  ZR_(az_, el_);

#undef ZR_HELPER_
#undef ZR_

}

bool Bearing::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoMeas.Bearing)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoff(127);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // optional double az = 1;
      case 1: {
        if (tag == 9) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &az_)));

        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(17)) goto parse_el;
        break;
      }

      // optional double el = 2;
      case 2: {
        if (tag == 17) {
         parse_el:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &el_)));

        } else {
          goto handle_unusual;
        }
        if (input->ExpectAtEnd()) goto success;
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0 ||
            ::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormatLite::SkipField(input, tag));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:ProtoMeas.Bearing)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:ProtoMeas.Bearing)
  return false;
#undef DO_
}

void Bearing::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:ProtoMeas.Bearing)
  // optional double az = 1;
  if (this->az() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(1, this->az(), output);
  }

  // optional double el = 2;
  if (this->el() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(2, this->el(), output);
  }

  // @@protoc_insertion_point(serialize_end:ProtoMeas.Bearing)
}

::google::protobuf::uint8* Bearing::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  // @@protoc_insertion_point(serialize_to_array_start:ProtoMeas.Bearing)
  // optional double az = 1;
  if (this->az() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(1, this->az(), target);
  }

  // optional double el = 2;
  if (this->el() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(2, this->el(), target);
  }

  // @@protoc_insertion_point(serialize_to_array_end:ProtoMeas.Bearing)
  return target;
}

int Bearing::ByteSize() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoMeas.Bearing)
  int total_size = 0;

  // optional double az = 1;
  if (this->az() != 0) {
    total_size += 1 + 8;
  }

  // optional double el = 2;
  if (this->el() != 0) {
    total_size += 1 + 8;
  }

  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void Bearing::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoMeas.Bearing)
  if (GOOGLE_PREDICT_FALSE(&from == this)) {
    ::google::protobuf::internal::MergeFromFail(__FILE__, __LINE__);
  }
  const Bearing* source = 
      ::google::protobuf::internal::DynamicCastToGenerated<const Bearing>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ProtoMeas.Bearing)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ProtoMeas.Bearing)
    MergeFrom(*source);
  }
}

void Bearing::MergeFrom(const Bearing& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ProtoMeas.Bearing)
  if (GOOGLE_PREDICT_FALSE(&from == this)) {
    ::google::protobuf::internal::MergeFromFail(__FILE__, __LINE__);
  }
  if (from.az() != 0) {
    set_az(from.az());
  }
  if (from.el() != 0) {
    set_el(from.el());
  }
}

void Bearing::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ProtoMeas.Bearing)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Bearing::CopyFrom(const Bearing& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ProtoMeas.Bearing)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Bearing::IsInitialized() const {

  return true;
}

void Bearing::Swap(Bearing* other) {
  if (other == this) return;
  InternalSwap(other);
}
void Bearing::InternalSwap(Bearing* other) {
  std::swap(az_, other->az_);
  std::swap(el_, other->el_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
  std::swap(_cached_size_, other->_cached_size_);
}

::google::protobuf::Metadata Bearing::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = Bearing_descriptor_;
  metadata.reflection = Bearing_reflection_;
  return metadata;
}

#if PROTOBUF_INLINE_NOT_IN_HEADERS
// Bearing

// optional double az = 1;
void Bearing::clear_az() {
  az_ = 0;
}
 double Bearing::az() const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Bearing.az)
  return az_;
}
 void Bearing::set_az(double value) {
  
  az_ = value;
  // @@protoc_insertion_point(field_set:ProtoMeas.Bearing.az)
}

// optional double el = 2;
void Bearing::clear_el() {
  el_ = 0;
}
 double Bearing::el() const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Bearing.el)
  return el_;
}
 void Bearing::set_el(double value) {
  
  el_ = value;
  // @@protoc_insertion_point(field_set:ProtoMeas.Bearing.el)
}

#endif  // PROTOBUF_INLINE_NOT_IN_HEADERS

// ===================================================================

#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Measurements::kNumFeaturePointsFieldNumber;
const int Measurements::kFeaturePointsFieldNumber;
const int Measurements::kBearingsFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Measurements::Measurements()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoMeas.Measurements)
}

void Measurements::InitAsDefaultInstance() {
  _is_default_instance_ = true;
}

Measurements::Measurements(const Measurements& from)
  : ::google::protobuf::Message(),
    _internal_metadata_(NULL) {
  SharedCtor();
  MergeFrom(from);
  // @@protoc_insertion_point(copy_constructor:ProtoMeas.Measurements)
}

void Measurements::SharedCtor() {
    _is_default_instance_ = false;
  _cached_size_ = 0;
  num_feature_points_ = 0u;
}

Measurements::~Measurements() {
  // @@protoc_insertion_point(destructor:ProtoMeas.Measurements)
  SharedDtor();
}

void Measurements::SharedDtor() {
  if (this != default_instance_) {
  }
}

void Measurements::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* Measurements::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return Measurements_descriptor_;
}

const Measurements& Measurements::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_measurement_2eproto();
  return *default_instance_;
}

Measurements* Measurements::default_instance_ = NULL;

Measurements* Measurements::New(::google::protobuf::Arena* arena) const {
  Measurements* n = new Measurements;
  if (arena != NULL) {
    arena->Own(n);
  }
  return n;
}

void Measurements::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoMeas.Measurements)
  num_feature_points_ = 0u;
  feature_points_.Clear();
  bearings_.Clear();
}

bool Measurements::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoMeas.Measurements)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoff(127);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // optional uint32 num_feature_points = 1;
      case 1: {
        if (tag == 8) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::uint32, ::google::protobuf::internal::WireFormatLite::TYPE_UINT32>(
                 input, &num_feature_points_)));

        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(18)) goto parse_feature_points;
        break;
      }

      // repeated .ProtoMeas.Position feature_points = 2;
      case 2: {
        if (tag == 18) {
         parse_feature_points:
          DO_(input->IncrementRecursionDepth());
         parse_loop_feature_points:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtualNoRecursionDepth(
                input, add_feature_points()));
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(18)) goto parse_loop_feature_points;
        if (input->ExpectTag(26)) goto parse_loop_bearings;
        input->UnsafeDecrementRecursionDepth();
        break;
      }

      // repeated .ProtoMeas.Bearing bearings = 3;
      case 3: {
        if (tag == 26) {
          DO_(input->IncrementRecursionDepth());
         parse_loop_bearings:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtualNoRecursionDepth(
                input, add_bearings()));
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(26)) goto parse_loop_bearings;
        input->UnsafeDecrementRecursionDepth();
        if (input->ExpectAtEnd()) goto success;
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0 ||
            ::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormatLite::SkipField(input, tag));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:ProtoMeas.Measurements)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:ProtoMeas.Measurements)
  return false;
#undef DO_
}

void Measurements::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:ProtoMeas.Measurements)
  // optional uint32 num_feature_points = 1;
  if (this->num_feature_points() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteUInt32(1, this->num_feature_points(), output);
  }

  // repeated .ProtoMeas.Position feature_points = 2;
  for (unsigned int i = 0, n = this->feature_points_size(); i < n; i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      2, this->feature_points(i), output);
  }

  // repeated .ProtoMeas.Bearing bearings = 3;
  for (unsigned int i = 0, n = this->bearings_size(); i < n; i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      3, this->bearings(i), output);
  }

  // @@protoc_insertion_point(serialize_end:ProtoMeas.Measurements)
}

::google::protobuf::uint8* Measurements::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  // @@protoc_insertion_point(serialize_to_array_start:ProtoMeas.Measurements)
  // optional uint32 num_feature_points = 1;
  if (this->num_feature_points() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteUInt32ToArray(1, this->num_feature_points(), target);
  }

  // repeated .ProtoMeas.Position feature_points = 2;
  for (unsigned int i = 0, n = this->feature_points_size(); i < n; i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageNoVirtualToArray(
        2, this->feature_points(i), false, target);
  }

  // repeated .ProtoMeas.Bearing bearings = 3;
  for (unsigned int i = 0, n = this->bearings_size(); i < n; i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageNoVirtualToArray(
        3, this->bearings(i), false, target);
  }

  // @@protoc_insertion_point(serialize_to_array_end:ProtoMeas.Measurements)
  return target;
}

int Measurements::ByteSize() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoMeas.Measurements)
  int total_size = 0;

  // optional uint32 num_feature_points = 1;
  if (this->num_feature_points() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::UInt32Size(
        this->num_feature_points());
  }

  // repeated .ProtoMeas.Position feature_points = 2;
  total_size += 1 * this->feature_points_size();
  for (int i = 0; i < this->feature_points_size(); i++) {
    total_size +=
      ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
        this->feature_points(i));
  }

  // repeated .ProtoMeas.Bearing bearings = 3;
  total_size += 1 * this->bearings_size();
  for (int i = 0; i < this->bearings_size(); i++) {
    total_size +=
      ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
        this->bearings(i));
  }

  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void Measurements::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoMeas.Measurements)
  if (GOOGLE_PREDICT_FALSE(&from == this)) {
    ::google::protobuf::internal::MergeFromFail(__FILE__, __LINE__);
  }
  const Measurements* source = 
      ::google::protobuf::internal::DynamicCastToGenerated<const Measurements>(
          &from);
  if (source == NULL) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ProtoMeas.Measurements)
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ProtoMeas.Measurements)
    MergeFrom(*source);
  }
}

void Measurements::MergeFrom(const Measurements& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ProtoMeas.Measurements)
  if (GOOGLE_PREDICT_FALSE(&from == this)) {
    ::google::protobuf::internal::MergeFromFail(__FILE__, __LINE__);
  }
  feature_points_.MergeFrom(from.feature_points_);
  bearings_.MergeFrom(from.bearings_);
  if (from.num_feature_points() != 0) {
    set_num_feature_points(from.num_feature_points());
  }
}

void Measurements::CopyFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ProtoMeas.Measurements)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Measurements::CopyFrom(const Measurements& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ProtoMeas.Measurements)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Measurements::IsInitialized() const {

  return true;
}

void Measurements::Swap(Measurements* other) {
  if (other == this) return;
  InternalSwap(other);
}
void Measurements::InternalSwap(Measurements* other) {
  std::swap(num_feature_points_, other->num_feature_points_);
  feature_points_.UnsafeArenaSwap(&other->feature_points_);
  bearings_.UnsafeArenaSwap(&other->bearings_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
  std::swap(_cached_size_, other->_cached_size_);
}

::google::protobuf::Metadata Measurements::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = Measurements_descriptor_;
  metadata.reflection = Measurements_reflection_;
  return metadata;
}

#if PROTOBUF_INLINE_NOT_IN_HEADERS
// Measurements

// optional uint32 num_feature_points = 1;
void Measurements::clear_num_feature_points() {
  num_feature_points_ = 0u;
}
 ::google::protobuf::uint32 Measurements::num_feature_points() const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Measurements.num_feature_points)
  return num_feature_points_;
}
 void Measurements::set_num_feature_points(::google::protobuf::uint32 value) {
  
  num_feature_points_ = value;
  // @@protoc_insertion_point(field_set:ProtoMeas.Measurements.num_feature_points)
}

// repeated .ProtoMeas.Position feature_points = 2;
int Measurements::feature_points_size() const {
  return feature_points_.size();
}
void Measurements::clear_feature_points() {
  feature_points_.Clear();
}
const ::ProtoMeas::Position& Measurements::feature_points(int index) const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Measurements.feature_points)
  return feature_points_.Get(index);
}
::ProtoMeas::Position* Measurements::mutable_feature_points(int index) {
  // @@protoc_insertion_point(field_mutable:ProtoMeas.Measurements.feature_points)
  return feature_points_.Mutable(index);
}
::ProtoMeas::Position* Measurements::add_feature_points() {
  // @@protoc_insertion_point(field_add:ProtoMeas.Measurements.feature_points)
  return feature_points_.Add();
}
::google::protobuf::RepeatedPtrField< ::ProtoMeas::Position >*
Measurements::mutable_feature_points() {
  // @@protoc_insertion_point(field_mutable_list:ProtoMeas.Measurements.feature_points)
  return &feature_points_;
}
const ::google::protobuf::RepeatedPtrField< ::ProtoMeas::Position >&
Measurements::feature_points() const {
  // @@protoc_insertion_point(field_list:ProtoMeas.Measurements.feature_points)
  return feature_points_;
}

// repeated .ProtoMeas.Bearing bearings = 3;
int Measurements::bearings_size() const {
  return bearings_.size();
}
void Measurements::clear_bearings() {
  bearings_.Clear();
}
const ::ProtoMeas::Bearing& Measurements::bearings(int index) const {
  // @@protoc_insertion_point(field_get:ProtoMeas.Measurements.bearings)
  return bearings_.Get(index);
}
::ProtoMeas::Bearing* Measurements::mutable_bearings(int index) {
  // @@protoc_insertion_point(field_mutable:ProtoMeas.Measurements.bearings)
  return bearings_.Mutable(index);
}
::ProtoMeas::Bearing* Measurements::add_bearings() {
  // @@protoc_insertion_point(field_add:ProtoMeas.Measurements.bearings)
  return bearings_.Add();
}
::google::protobuf::RepeatedPtrField< ::ProtoMeas::Bearing >*
Measurements::mutable_bearings() {
  // @@protoc_insertion_point(field_mutable_list:ProtoMeas.Measurements.bearings)
  return &bearings_;
}
const ::google::protobuf::RepeatedPtrField< ::ProtoMeas::Bearing >&
Measurements::bearings() const {
  // @@protoc_insertion_point(field_list:ProtoMeas.Measurements.bearings)
  return bearings_;
}

#endif  // PROTOBUF_INLINE_NOT_IN_HEADERS

// @@protoc_insertion_point(namespace_scope)

}  // namespace ProtoMeas

// @@protoc_insertion_point(global_scope)