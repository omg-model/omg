function [ diag_fcns ] = diagnostic_functions ( )
    
    diag_fcns.f2v=@f2v;
    diag_fcns.v2f=@v2f;
    diag_fcns.read_genie_netcdf=@read_genie_netcdf;
    diag_fcns.read_genie_netcdf_slice=@read_genie_netcdf_slice;
    diag_fcns.generate_save_indices=@generate_save_indices;
    diag_fcns.integrate_output=@integrate_output;
    diag_fcns.create_netcdf_output=@create_netcdf_output;
    diag_fcns.write_netcdf_output=@write_netcdf_output;
    diag_fcns.write_sparse_output=@write_sparse_output;
    diag_fcns.reset_output=@reset_output;
    diag_fcns.mol_clock=@mol_clock;
    diag_fcns.gene_mutate=@gene_mutate;
    diag_fcns.maxrows=@maxrows;
    diag_fcns.maxcols=@maxcols;

end


function [ data ] = load_variable(matObj,variable,nyr)
% extract output data from matObj file
% variables:
%     {'AGE'     }
%     {'DOP'     }
%     {'GENOME'  }
%     {'PHY'     }
%     {'PO4'     }

data = cell2mat(matObj.(variable)(nyr,1));

