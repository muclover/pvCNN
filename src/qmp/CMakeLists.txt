include_directories(.)

add_library(qmp
  commitment.cpp
  matrix_gadget.cpp
  model.cpp
  timer.cpp
  utils.cpp
  writer.cpp
)

target_link_libraries(
  qmp

  snark
  cryptopp
)

target_include_directories(
  qmp

  PUBLIC
  ${DEPENDS_DIR}/libsnark
  ${DEPENDS_DIR}/libsnark/depends/libfqfft
)