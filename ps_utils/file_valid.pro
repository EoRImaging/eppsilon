function file_valid, filename
  return, file_test(filename) *  (1 - file_test(filename, /zero_length))
end
