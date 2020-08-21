// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: measurement.proto

#include "measurement.pb.h"

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

namespace protobuf_measurement_2eproto {
extern PROTOBUF_INTERNAL_EXPORT_protobuf_measurement_2eproto ::google::protobuf::internal::SCCInfo<0> scc_info_Bearing;
extern PROTOBUF_INTERNAL_EXPORT_protobuf_measurement_2eproto ::google::protobuf::internal::SCCInfo<0> scc_info_Position;
}  // namespace protobuf_measurement_2eproto
namespace ProtoMeas {
class PositionDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<Position>
      _instance;
} _Position_default_instance_;
class BearingDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<Bearing>
      _instance;
} _Bearing_default_instance_;
class MeasurementsDefaultTypeInternal {
 public:
  ::google::protobuf::internal::ExplicitlyConstructed<Measurements>
      _instance;
} _Measurements_default_instance_;
}  // namespace ProtoMeas
namespace protobuf_measurement_2eproto {
static void InitDefaultsPosition() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ProtoMeas::_Position_default_instance_;
    new (ptr) ::ProtoMeas::Position();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ProtoMeas::Position::InitAsDefaultInstance();
}

::google::protobuf::internal::SCCInfo<0> scc_info_Position =
    {{ATOMIC_VAR_INIT(::google::protobuf::internal::SCCInfoBase::kUninitialized), 0, InitDefaultsPosition}, {}};

static void InitDefaultsBearing() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ProtoMeas::_Bearing_default_instance_;
    new (ptr) ::ProtoMeas::Bearing();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ProtoMeas::Bearing::InitAsDefaultInstance();
}

::google::protobuf::internal::SCCInfo<0> scc_info_Bearing =
    {{ATOMIC_VAR_INIT(::google::protobuf::internal::SCCInfoBase::kUninitialized), 0, InitDefaultsBearing}, {}};

static void InitDefaultsMeasurements() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ProtoMeas::_Measurements_default_instance_;
    new (ptr) ::ProtoMeas::Measurements();
    ::google::protobuf::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ProtoMeas::Measurements::InitAsDefaultInstance();
}

::google::protobuf::internal::SCCInfo<2> scc_info_Measurements =
    {{ATOMIC_VAR_INIT(::google::protobuf::internal::SCCInfoBase::kUninitialized), 2, InitDefaultsMeasurements}, {
      &protobuf_measurement_2eproto::scc_info_Position.base,
      &protobuf_measurement_2eproto::scc_info_Bearing.base,}};

void InitDefaults() {
  ::google::protobuf::internal::InitSCC(&scc_info_Position.base);
  ::google::protobuf::internal::InitSCC(&scc_info_Bearing.base);
  ::google::protobuf::internal::InitSCC(&scc_info_Measurements.base);
}

::google::protobuf::Metadata file_level_metadata[3];

const ::google::protobuf::uint32 TableStruct::offsets[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Position, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Position, x_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Position, y_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Position, z_),
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Bearing, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Bearing, az_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Bearing, el_),
  ~0u,  // no _has_bits_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Measurements, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Measurements, num_feature_points_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Measurements, feature_points_),
  GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(::ProtoMeas::Measurements, bearings_),
};
static const ::google::protobuf::internal::MigrationSchema schemas[] GOOGLE_PROTOBUF_ATTRIBUTE_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, sizeof(::ProtoMeas::Position)},
  { 8, -1, sizeof(::ProtoMeas::Bearing)},
  { 15, -1, sizeof(::ProtoMeas::Measurements)},
};

static ::google::protobuf::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::google::protobuf::Message*>(&::ProtoMeas::_Position_default_instance_),
  reinterpret_cast<const ::google::protobuf::Message*>(&::ProtoMeas::_Bearing_default_instance_),
  reinterpret_cast<const ::google::protobuf::Message*>(&::ProtoMeas::_Measurements_default_instance_),
};

void protobuf_AssignDescriptors() {
  AddDescriptors();
  AssignDescriptors(
      "measurement.proto", schemas, file_default_instances, TableStruct::offsets,
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
      "\n\021measurement.proto\022\tProtoMeas\"+\n\010Positi"
      "on\022\t\n\001x\030\001 \001(\001\022\t\n\001y\030\002 \001(\001\022\t\n\001z\030\003 \001(\001\"!\n\007B"
      "earing\022\n\n\002az\030\001 \001(\001\022\n\n\002el\030\002 \001(\001\"}\n\014Measur"
      "ements\022\032\n\022num_feature_points\030\001 \001(\r\022+\n\016fe"
      "ature_points\030\002 \003(\0132\023.ProtoMeas.Position\022"
      "$\n\010bearings\030\003 \003(\0132\022.ProtoMeas.Bearingb\006p"
      "roto3"
  };
  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
      descriptor, 245);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "measurement.proto", &protobuf_RegisterTypes);
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
}  // namespace protobuf_measurement_2eproto
namespace ProtoMeas {

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
      &protobuf_measurement_2eproto::scc_info_Position.base);
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoMeas.Position)
}
Position::Position(const Position& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::memcpy(&x_, &from.x_,
    static_cast<size_t>(reinterpret_cast<char*>(&z_) -
    reinterpret_cast<char*>(&x_)) + sizeof(z_));
  // @@protoc_insertion_point(copy_constructor:ProtoMeas.Position)
}

void Position::SharedCtor() {
  ::memset(&x_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&z_) -
      reinterpret_cast<char*>(&x_)) + sizeof(z_));
}

Position::~Position() {
  // @@protoc_insertion_point(destructor:ProtoMeas.Position)
  SharedDtor();
}

void Position::SharedDtor() {
}

void Position::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ::google::protobuf::Descriptor* Position::descriptor() {
  ::protobuf_measurement_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_measurement_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const Position& Position::default_instance() {
  ::google::protobuf::internal::InitSCC(&protobuf_measurement_2eproto::scc_info_Position.base);
  return *internal_default_instance();
}


void Position::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoMeas.Position)
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
  // @@protoc_insertion_point(parse_start:ProtoMeas.Position)
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
  // @@protoc_insertion_point(serialize_end:ProtoMeas.Position)
}

::google::protobuf::uint8* Position::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:ProtoMeas.Position)
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
  // @@protoc_insertion_point(serialize_to_array_end:ProtoMeas.Position)
  return target;
}

size_t Position::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoMeas.Position)
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
// @@protoc_insertion_point(generalized_merge_from_start:ProtoMeas.Position)
  GOOGLE_DCHECK_NE(&from, this);
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
  using std::swap;
  swap(x_, other->x_);
  swap(y_, other->y_);
  swap(z_, other->z_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
}

::google::protobuf::Metadata Position::GetMetadata() const {
  protobuf_measurement_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_measurement_2eproto::file_level_metadata[kIndexInFileMessages];
}


// ===================================================================

void Bearing::InitAsDefaultInstance() {
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Bearing::kAzFieldNumber;
const int Bearing::kElFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Bearing::Bearing()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  ::google::protobuf::internal::InitSCC(
      &protobuf_measurement_2eproto::scc_info_Bearing.base);
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoMeas.Bearing)
}
Bearing::Bearing(const Bearing& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::memcpy(&az_, &from.az_,
    static_cast<size_t>(reinterpret_cast<char*>(&el_) -
    reinterpret_cast<char*>(&az_)) + sizeof(el_));
  // @@protoc_insertion_point(copy_constructor:ProtoMeas.Bearing)
}

void Bearing::SharedCtor() {
  ::memset(&az_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&el_) -
      reinterpret_cast<char*>(&az_)) + sizeof(el_));
}

Bearing::~Bearing() {
  // @@protoc_insertion_point(destructor:ProtoMeas.Bearing)
  SharedDtor();
}

void Bearing::SharedDtor() {
}

void Bearing::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ::google::protobuf::Descriptor* Bearing::descriptor() {
  ::protobuf_measurement_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_measurement_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const Bearing& Bearing::default_instance() {
  ::google::protobuf::internal::InitSCC(&protobuf_measurement_2eproto::scc_info_Bearing.base);
  return *internal_default_instance();
}


void Bearing::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoMeas.Bearing)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  ::memset(&az_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&el_) -
      reinterpret_cast<char*>(&az_)) + sizeof(el_));
  _internal_metadata_.Clear();
}

bool Bearing::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoMeas.Bearing)
  for (;;) {
    ::std::pair<::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // double az = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(9u /* 9 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &az_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // double el = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(17u /* 17 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, &el_)));
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
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double az = 1;
  if (this->az() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(1, this->az(), output);
  }

  // double el = 2;
  if (this->el() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(2, this->el(), output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:ProtoMeas.Bearing)
}

::google::protobuf::uint8* Bearing::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:ProtoMeas.Bearing)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // double az = 1;
  if (this->az() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(1, this->az(), target);
  }

  // double el = 2;
  if (this->el() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteDoubleToArray(2, this->el(), target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ProtoMeas.Bearing)
  return target;
}

size_t Bearing::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoMeas.Bearing)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // double az = 1;
  if (this->az() != 0) {
    total_size += 1 + 8;
  }

  // double el = 2;
  if (this->el() != 0) {
    total_size += 1 + 8;
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void Bearing::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoMeas.Bearing)
  GOOGLE_DCHECK_NE(&from, this);
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
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

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
  using std::swap;
  swap(az_, other->az_);
  swap(el_, other->el_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
}

::google::protobuf::Metadata Bearing::GetMetadata() const {
  protobuf_measurement_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_measurement_2eproto::file_level_metadata[kIndexInFileMessages];
}


// ===================================================================

void Measurements::InitAsDefaultInstance() {
}
#if !defined(_MSC_VER) || _MSC_VER >= 1900
const int Measurements::kNumFeaturePointsFieldNumber;
const int Measurements::kFeaturePointsFieldNumber;
const int Measurements::kBearingsFieldNumber;
#endif  // !defined(_MSC_VER) || _MSC_VER >= 1900

Measurements::Measurements()
  : ::google::protobuf::Message(), _internal_metadata_(NULL) {
  ::google::protobuf::internal::InitSCC(
      &protobuf_measurement_2eproto::scc_info_Measurements.base);
  SharedCtor();
  // @@protoc_insertion_point(constructor:ProtoMeas.Measurements)
}
Measurements::Measurements(const Measurements& from)
  : ::google::protobuf::Message(),
      _internal_metadata_(NULL),
      feature_points_(from.feature_points_),
      bearings_(from.bearings_) {
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  num_feature_points_ = from.num_feature_points_;
  // @@protoc_insertion_point(copy_constructor:ProtoMeas.Measurements)
}

void Measurements::SharedCtor() {
  num_feature_points_ = 0u;
}

Measurements::~Measurements() {
  // @@protoc_insertion_point(destructor:ProtoMeas.Measurements)
  SharedDtor();
}

void Measurements::SharedDtor() {
}

void Measurements::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ::google::protobuf::Descriptor* Measurements::descriptor() {
  ::protobuf_measurement_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_measurement_2eproto::file_level_metadata[kIndexInFileMessages].descriptor;
}

const Measurements& Measurements::default_instance() {
  ::google::protobuf::internal::InitSCC(&protobuf_measurement_2eproto::scc_info_Measurements.base);
  return *internal_default_instance();
}


void Measurements::Clear() {
// @@protoc_insertion_point(message_clear_start:ProtoMeas.Measurements)
  ::google::protobuf::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  feature_points_.Clear();
  bearings_.Clear();
  num_feature_points_ = 0u;
  _internal_metadata_.Clear();
}

bool Measurements::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!GOOGLE_PREDICT_TRUE(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:ProtoMeas.Measurements)
  for (;;) {
    ::std::pair<::google::protobuf::uint32, bool> p = input->ReadTagWithCutoffNoLastTag(127u);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // uint32 num_feature_points = 1;
      case 1: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(8u /* 8 & 0xFF */)) {

          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::uint32, ::google::protobuf::internal::WireFormatLite::TYPE_UINT32>(
                 input, &num_feature_points_)));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // repeated .ProtoMeas.Position feature_points = 2;
      case 2: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(18u /* 18 & 0xFF */)) {
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessage(
                input, add_feature_points()));
        } else {
          goto handle_unusual;
        }
        break;
      }

      // repeated .ProtoMeas.Bearing bearings = 3;
      case 3: {
        if (static_cast< ::google::protobuf::uint8>(tag) ==
            static_cast< ::google::protobuf::uint8>(26u /* 26 & 0xFF */)) {
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessage(
                input, add_bearings()));
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
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // uint32 num_feature_points = 1;
  if (this->num_feature_points() != 0) {
    ::google::protobuf::internal::WireFormatLite::WriteUInt32(1, this->num_feature_points(), output);
  }

  // repeated .ProtoMeas.Position feature_points = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->feature_points_size()); i < n; i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      2,
      this->feature_points(static_cast<int>(i)),
      output);
  }

  // repeated .ProtoMeas.Bearing bearings = 3;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->bearings_size()); i < n; i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      3,
      this->bearings(static_cast<int>(i)),
      output);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), output);
  }
  // @@protoc_insertion_point(serialize_end:ProtoMeas.Measurements)
}

::google::protobuf::uint8* Measurements::InternalSerializeWithCachedSizesToArray(
    bool deterministic, ::google::protobuf::uint8* target) const {
  (void)deterministic; // Unused
  // @@protoc_insertion_point(serialize_to_array_start:ProtoMeas.Measurements)
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // uint32 num_feature_points = 1;
  if (this->num_feature_points() != 0) {
    target = ::google::protobuf::internal::WireFormatLite::WriteUInt32ToArray(1, this->num_feature_points(), target);
  }

  // repeated .ProtoMeas.Position feature_points = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->feature_points_size()); i < n; i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageToArray(
        2, this->feature_points(static_cast<int>(i)), deterministic, target);
  }

  // repeated .ProtoMeas.Bearing bearings = 3;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->bearings_size()); i < n; i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      InternalWriteMessageToArray(
        3, this->bearings(static_cast<int>(i)), deterministic, target);
  }

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ProtoMeas.Measurements)
  return target;
}

size_t Measurements::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ProtoMeas.Measurements)
  size_t total_size = 0;

  if ((_internal_metadata_.have_unknown_fields() &&  ::google::protobuf::internal::GetProto3PreserveUnknownsDefault())) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        (::google::protobuf::internal::GetProto3PreserveUnknownsDefault()   ? _internal_metadata_.unknown_fields()   : _internal_metadata_.default_instance()));
  }
  // repeated .ProtoMeas.Position feature_points = 2;
  {
    unsigned int count = static_cast<unsigned int>(this->feature_points_size());
    total_size += 1UL * count;
    for (unsigned int i = 0; i < count; i++) {
      total_size +=
        ::google::protobuf::internal::WireFormatLite::MessageSize(
          this->feature_points(static_cast<int>(i)));
    }
  }

  // repeated .ProtoMeas.Bearing bearings = 3;
  {
    unsigned int count = static_cast<unsigned int>(this->bearings_size());
    total_size += 1UL * count;
    for (unsigned int i = 0; i < count; i++) {
      total_size +=
        ::google::protobuf::internal::WireFormatLite::MessageSize(
          this->bearings(static_cast<int>(i)));
    }
  }

  // uint32 num_feature_points = 1;
  if (this->num_feature_points() != 0) {
    total_size += 1 +
      ::google::protobuf::internal::WireFormatLite::UInt32Size(
        this->num_feature_points());
  }

  int cached_size = ::google::protobuf::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void Measurements::MergeFrom(const ::google::protobuf::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ProtoMeas.Measurements)
  GOOGLE_DCHECK_NE(&from, this);
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
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom(from._internal_metadata_);
  ::google::protobuf::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

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
  using std::swap;
  CastToBase(&feature_points_)->InternalSwap(CastToBase(&other->feature_points_));
  CastToBase(&bearings_)->InternalSwap(CastToBase(&other->bearings_));
  swap(num_feature_points_, other->num_feature_points_);
  _internal_metadata_.Swap(&other->_internal_metadata_);
}

::google::protobuf::Metadata Measurements::GetMetadata() const {
  protobuf_measurement_2eproto::protobuf_AssignDescriptorsOnce();
  return ::protobuf_measurement_2eproto::file_level_metadata[kIndexInFileMessages];
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace ProtoMeas
namespace google {
namespace protobuf {
template<> GOOGLE_PROTOBUF_ATTRIBUTE_NOINLINE ::ProtoMeas::Position* Arena::CreateMaybeMessage< ::ProtoMeas::Position >(Arena* arena) {
  return Arena::CreateInternal< ::ProtoMeas::Position >(arena);
}
template<> GOOGLE_PROTOBUF_ATTRIBUTE_NOINLINE ::ProtoMeas::Bearing* Arena::CreateMaybeMessage< ::ProtoMeas::Bearing >(Arena* arena) {
  return Arena::CreateInternal< ::ProtoMeas::Bearing >(arena);
}
template<> GOOGLE_PROTOBUF_ATTRIBUTE_NOINLINE ::ProtoMeas::Measurements* Arena::CreateMaybeMessage< ::ProtoMeas::Measurements >(Arena* arena) {
  return Arena::CreateInternal< ::ProtoMeas::Measurements >(arena);
}
}  // namespace protobuf
}  // namespace google

// @@protoc_insertion_point(global_scope)
