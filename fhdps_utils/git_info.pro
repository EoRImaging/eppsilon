function git_info, repo_path

  git, 'config --get remote.origin.url', repo_path=repo_path, result=origin
  git, 'rev-parse HEAD', repo_path=repo_path, result=hash
  git, 'describe --dirty --tag --always', repo_path=repo_path, result=description
  git, 'rev-parse --abbrev-ref HEAD', repo_path=repo_path, result=branch

  ;; turn them into strings instead of arrays
  origin = origin[0]
  branch = branch[0]
  description = description[0]
  hash = hash[0]

  git_string = origin + '_' + branch + '_' + description + '_' + hash

  return, {origin: origin, branch: branch, description: description, hash: hash, $
           full_string: git_string}
end
