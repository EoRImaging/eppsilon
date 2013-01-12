function getvar_savefile, savefile, varname
  savefile_obj = obj_new('idl_savefile', savefile)
  savefile_obj->Restore, varname
  obj_destroy, savefile_obj

  q=execute('return,'+varname) 
end

