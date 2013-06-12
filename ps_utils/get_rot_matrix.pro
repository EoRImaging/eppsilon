function get_rot_matrix, theta, phi, inverse = inverse

  rot_matrix = [[cos(theta) * cos(phi)^2d + sin(phi)^2d, (cos(theta)-1d) * cos(phi)*sin(phi), sin(theta) * cos(phi)], $
                [(cos(theta)-1d) * cos(phi)*sin(phi), cos(theta) * sin(phi)^2d + cos(phi)^2d, sin(theta) * sin(phi)], $
                [(-1) * sin(theta) * cos(phi), (-1) * sin(theta) * sin(phi), cos(theta)]]

  if keyword_set(inverse) then return, invert(rot_matrix) else return, rot_matrix

end
