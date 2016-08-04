function out = casadi2yalmip(nlp,data,discrete)
  import casadi.*
  fun = Function('nlp',nlp,char('x','p'),char('f','g'));
  fun = fun.expand();

  ins = fun.sx_in();
  fun = Function('f',ins,fun.call(ins),struct('live_variables',false));

  lbg = -inf*DM(fun.sparsity_out(1));
  ubg = inf*DM(fun.sparsity_out(1));
  if isfield(data,'lbg'), lbg(:,:) = data.lbg; end
  if isfield(data,'ubg'), ubg(:,:) = data.ubg; end


  lbx = -inf*DM(fun.sparsity_in(0));
  x0 = DM(fun.sparsity_in(0));
  ubx = inf*DM(fun.sparsity_in(0));
  if isfield(data,'lbx'), lbx(:,:) = data.lbx; end
  if isfield(data,'ubx'), ubx(:,:) = data.ubx; end
  if isfield(data,'x0'), x0(:,:) = data.x0; end
  
  p = DM(fun.sparsity_in(1));
  if isfield(data,'p'), p(:,:) = data.p; end
  
  main = fopen('exported.m','w');
  
  for i=1:size(x0,1)
    if discrete(i)
        fprintf(main,'x%d = intvar(1,1);assign(x%d,%.16f);\n',i-1,i-1,full(x0(i)));
    else
        fprintf(main,'x%d = sdpvar(1,1);assign(x%d,%.16f);\n',i-1,i-1,full(x0(i)));
    end
  end

  for i=1:size(p,1)
    fprintf(main,'p%d = %.16f;\n',i-1,full(p(i)));
  end

  algorithm = strsplit(fun.getDescription(),'\n');
  for i=1:numel(algorithm)
    l = algorithm{i};
    if strfind(l, '@')
      break
    end
  end
  algorithm = strjoin(algorithm(i:end),'\n');
    
  algorithm = regexprep(algorithm,'@','at');
  algorithm = regexprep(algorithm,'input\[0\]\[(\d+)\]','x$1');
  algorithm = regexprep(algorithm,'input\[1\]\[(\d+)\]','p$1');
  algorithm = regexprep(algorithm,'output\[0\]\[(\d+)\]','f');
  algorithm = regexprep(algorithm,'output\[1\]\[(\d+)\]','g$1');
  algorithm = regexprep(algorithm,'sq\((.*?)\)','(\1)^2');

  fprintf(main,algorithm);
  
  fprintf(main,'g = [');

  for i=1:size(lbx,1)
    
    if isinf(full(lbx(i))) && isinf(full(ubx(i)))
      continue;
    elseif full(lbx(i))==full(ubx(i))
      fprintf(main,'x%i == %.16f',i-1,full(lbx(i)));
    elseif isinf(full(lbx(i)))
      fprintf(main,'x%i <= %.16f',i-1,full(ubx(i)));
    elseif isinf(full(ubx(i)))
      fprintf(main,'x%i >= %.16f',i-1,full(lbx(i)));
    else
      fprintf(main,'x%i >= %.16f,',i-1,full(lbx(i)));
      fprintf(main,'x%i <= %.16f',i-1,full(ubx(i)));
    end
    fprintf(main,',...\n');
  end
  
  for i=1:size(lbg,1)
    if isinf(full(lbg(i))) && isinf(full(ubg(i)))
      continue;
    elseif full(lbg(i))==full(ubg(i))
      fprintf(main,'g%i == %.16f',i-1,full(lbg(i)));
    elseif isinf(full(lbg(i)))
      fprintf(main,'g%i <= %.16f',i-1,full(ubg(i)));
    elseif isinf(full(ubg(i)))
      fprintf(main,'g%i >= %.16f',i-1,full(lbg(i)));
    else
      fprintf(main,'g%i >= %.16f,',i-1,full(lbg(i)));
      fprintf(main,'g%i <= %.16f',i-1,full(ubg(i)));
    end
    fprintf(main,',...\n');
  end
  
  fprintf(main,'];\n');
  fprintf(main,'ops = sdpsettings(''solver'',''bonmin'');');
  fprintf(main,'optimize(g,f,ops)');
end
