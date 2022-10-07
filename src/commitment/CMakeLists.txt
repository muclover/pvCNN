include_directories(.)

find_package(SEAL 4.0 REQUIRED)

add_executable(
  HE-main

  he-main.cpp
)

target_link_libraries(
  HE-main

  snark
  SEAL::seal
)

target_include_directories(
  HE-main

  PUBLIC
  ${DEPENDS_DIR}/libsnark
  ${DEPENDS_DIR}/libsnark/depends/libfqfft
)
