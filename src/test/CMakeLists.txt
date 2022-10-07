include_directories(.)

add_executable(
  test_matrix_MNIST

  matrix.cpp
)
target_link_libraries(
  test_matrix_MNIST

  qmp
)

add_executable(
  QAP_Matrix

  matrix_only_QAP.cpp
)
target_link_libraries(
  QAP_Matrix

  qmp
)

add_executable(
  QMP_Matrix

  matrix_only_QMP.cpp
)
target_link_libraries(
  QMP_Matrix

  qmp
)
add_executable(
  Non_linear

  test_non_linear_operation.cpp
)
target_link_libraries(
  Non_linear
  qmp
)
add_executable(
  model_lenet

  model_lenet.cpp
)
target_link_libraries(
  model_lenet

  qmp
)