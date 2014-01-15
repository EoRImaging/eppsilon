function getvar_savefile, savefile, varname, pointer_return=pointer_return,verbose=verbose
  savefile_obj = obj_new('idl_savefile', savefile)
  savefile_obj->Restore, varname
  obj_destroy, savefile_obj

  IF Keyword_Set(verbose) THEN print,"Restoring "+varname+" from file: "+file_basename(savefile,'.sav',/fold_case)
  IF Keyword_Set(pointer_return) THEN p=execute(varname+'=Ptr_new('+varname+')')
  q=execute('return,'+varname) 
end

