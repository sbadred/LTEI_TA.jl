#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

    """These functions read files in xyz form and convert them to useful data
    for our test cases"""

    import XLSX
    struct molecule
        molecule_name::String
        Number_atoms::Int64
        basis_set_name::String
        Nb::Int64
    end
    struct orbitals
        name :: String
        expo
        coeff
        couche
    end
    struct atom
        Geo
        Orbits:: Array{orbitals, 1}
    end

    function mol_prop(file_name::String)
        xf = XLSX.readxlsx(file_name)
        sheetnames=XLSX.sheetnames(xf);
        #extract first sheet
        x1=xf[1]
        M=XLSX.getdata(x1);
        #Defining molecule structure
        Mol=molecule(x1[1,1],x1[1,2],x1[2,1],x1[3,1]);
        arrayOfAtoms = Vector{atom}(undef, x1[1,2])
        k=0;
        for i in 2:length(sheetnames)
            x2=xf[i];
            M1=XLSX.getdata(x2);
            natoms=x2[1,2];
            for j in 1:natoms
                k=k+1;
                valuetobefound=string(x2[1,1],string(j))
                Index=findfirst((x -> x==valuetobefound), M);
                Geom=[parse(Float64, x1[Index[1],Index[2]+1]) parse(Float64, x1[Index[1],Index[2]+2]) parse(Float64, x1[Index[1],Index[2]+3]) ]
                #Get each atom's orbitals info
                m=1;l=3;
                arrayOforbitals = Vector{orbitals}(undef, x2[1,3])
                while m <= x2[1,3]
                    expon=parse.(Float64,M1[l:l+M1[l-1,2]-1,2]);
                    coeff=parse.(Float64,M1[l:l+x2[l-1,2]-1,3]);

                    #Orbital S
                    if cmp(x2[l-1,1],"S")==0
                        arrayOforbitals[m]=orbitals("S",expon,(((2*expon)/pi).^(3/4)).*coeff,[0 0 0])
                        m=m+1;
                    end

                    #Orbital P
                    if cmp(x2[l-1,1],"P")==0
                        #Px consrtuct
                        arrayOforbitals[m]=orbitals("Px",expon,(((128*(expon.^5))./(pi^3)).^(1/4)).*coeff,[1 0 0])
                        m=m+1;
                        #Py construct
                        arrayOforbitals[m]=orbitals("Py",expon,(((128*(expon.^5))./(pi^3)).^(1/4)).*coeff,[0 1 0])
                        m=m+1;
                        #Pz construct
                        arrayOforbitals[m]=orbitals("Pz",expon,(((128*(expon.^5))./(pi^3)).^(1/4)).*coeff,[0 0 1])
                        m=m+1;
                    end

                    #Orbital D
                    if cmp(x2[l-1,1],"D")==0
                        #Dxx construxt
                        arrayOforbitals[m]=orbitals("Dxx",expon,(((2048*(expon.^7))./(9*pi^3)).^(1/4)).*coeff,[2 0 0])
                        m=m+1
                        #Dxy construct
                        arrayOforbitals[m]=orbitals("Dxy",expon,(((2048*(expon.^7))./(pi^3)).^(1/4)).*coeff,[1 1 0])
                        m=m+1
                        #Dxz construct
                        arrayOforbitals[m]=orbitals("Dxz",expon,(((2048*(expon.^7))./(pi^3)).^(1/4)).*coeff,[1 0 1])
                        m=m+1
                        #Dyy construct
                        arrayOforbitals[m]=orbitals("Dyy",expon,(((2048*(expon.^7))./(9*pi^3)).^(1/4)).*coeff,[0 2 0])
                        m=m+1
                        #Dyz construct
                        arrayOforbitals[m]=orbitals("Dyz",expon,(((2048*(expon.^7))./(pi^3)).^(1/4)).*coeff,[0 1 1])
                        m=m+1
                        #Dzz
                        arrayOforbitals[m]=orbitals("Dxy",expon,(((2048*(expon.^7))./(9*pi^3)).^(1/4)).*coeff,[0 0 2])
                        m=m+1
                    end

                    #Orbital F
                    if cmp(x2[l-1,1],"F")==0
                        #Fx construct
                        arrayOforbitals[m]=orbitals("Fxx",expon,(((32768*(expon.^9))./(225*pi^3)).^(1/4)).*coeff,[3 0 0])
                        m=m+1
                        #Fx2y
                        arrayOforbitals[m]=orbitals("Fx2y",expon,(((32768*(expon.^9))./(9*pi^3)).^(1/4)).*coeff,[2 1 0])
                        m=m+1
                        #Fx2z
                        arrayOforbitals[m]=orbitals("Fx2z",expon,(((32768*(expon.^9))./(9*pi^3)).^(1/4)).*coeff,[2 0 1])
                        m=m+1
                        #Fxy2
                        arrayOforbitals[m]=orbitals("Fxy2",expon,(((32768*(expon.^9))./(9*pi^3)).^(1/4)).*coeff,[1 2 0])
                        m=m+1
                        #Fxyz
                        arrayOforbitals[m]=orbitals("Fxyz",expon,(((32768*(expon.^9))./(pi^3)).^(1/4)).*coeff,[1 1 1])
                        m=m+1
                        #Fxz2
                        arrayOforbitals[m]=orbitals("Fxz2",expon,(((32768*(expon.^9))./(9*pi^3)).^(1/4)).*coeff,[1 0 2])
                        m=m+1
                        #Fy3
                        arrayOforbitals[m]=orbitals("Fy3",expon,(((32768*(expon.^9))./(225*pi^3)).^(1/4)).*coeff,[0 3 0])
                        m=m+1
                        #Fy2z
                        arrayOforbitals[m]=orbitals("Fxy2z",expon,(((32768*(expon.^9))./(9*pi^3)).^(1/4)).*coeff,[0 2 1])
                        m=m+1
                        #Fyz2
                        arrayOforbitals[m]=orbitals("Fyz2",expon,(((32768*(expon.^9))./(9*pi^3)).^(1/4)).*coeff,[0 1 2])
                        m=m+1
                        #Fz3
                        arrayOforbitals[m]=orbitals("Fz3",expon,(((32768*(expon.^9))./(225*pi^3)).^(1/4)).*coeff,[0 0 3])
                        m=m+1
                    end

                    #G orbitals
                    if cmp(x2[l-1,1],"G")==0
                        #Gx4
                        arrayOforbitals[m]=orbitals("Gx4",expon,(((524288*(expon.^11))./(11025*pi^3)).^(1/4)).*coeff,[4 0 0])
                        m=m+1
                        #Gx3y
                        arrayOforbitals[m]=orbitals("Gx3y",expon,(((524288*(expon.^11))./(225*pi^3)).^(1/4)).*coeff,[3 1 0])
                        m=m+1
                        #Gx3z
                        arrayOforbitals[m]=orbitals("Gx3z",expon,(((524288*(expon.^11))./(225*pi^3)).^(1/4)).*coeff,[3 0 1])
                        m=m+1
                        #Gx2y2
                        arrayOforbitals[m]=orbitals("Gx2y2",expon,(((524288*(expon.^11))./(81*pi^3)).^(1/4)).*coeff,[2 2 0])
                        m=m+1
                        #Gx2yz
                        arrayOforbitals[m]=orbitals("Gx2yz",expon,(((524288*(expon.^11))./(9*pi^3)).^(1/4)).*coeff,[2 1 1])
                        m=m+1
                        #Gx2z2
                        arrayOforbitals[m]=orbitals("Gx2z2",expon,(((524288*(expon.^11))./(81*pi^3)).^(1/4)).*coeff,[2 0 2])
                        m=m+1
                        #Gxy3
                        arrayOforbitals[m]=orbitals("Gxy3",expon,(((524288*(expon.^11))./(225*pi^3)).^(1/4)).*coeff,[1 3 0])
                        m=m+1
                        #Gxy2z
                        arrayOforbitals[m]=orbitals("Gxy2z",expon,(((524288*(expon.^11))./(9*pi^3)).^(1/4)).*coeff,[1 2 1])
                        m=m+1
                        #Gxyz2
                        arrayOforbitals[m]=orbitals("Gxyz2",expon,(((524288*(expon.^11))./(9*pi^3)).^(1/4)).*coeff,[1 1 2])
                        m=m+1
                        #Gxz3
                        arrayOforbitals[m]=orbitals("Gxz3",expon,(((524288*(expon.^11))./(225*pi^3)).^(1/4)).*coeff,[1 0 3])
                        m=m+1
                        #Gy4
                        arrayOforbitals[m]=orbitals("Gx4",expon,(((524288*(expon.^11))./(11025*pi^3)).^(1/4)).*coeff,[0 4 0])
                        m=m+1
                        #Gy3z
                        arrayOforbitals[m]=orbitals("Gy3z",expon,(((524288*(expon.^11))./(225*pi^3)).^(1/4)).*coeff,[0 3 1])
                        m=m+1
                        #Gy2z2
                        arrayOforbitals[m]=orbitals("Gy2z2",expon,(((524288*(expon.^11))./(81*pi^3)).^(1/4)).*coeff,[0 2 2])
                        m=m+1
                        #Gyz3
                        arrayOforbitals[m]=orbitals("Gyz3",expon,(((524288*(expon.^11))./(225*pi^3)).^(1/4)).*coeff,[0 1 3])
                        m=m+1
                        #Gz4
                        arrayOforbitals[m]=orbitals("Gz4",expon,(((524288*(expon.^11))./(11025*pi^3)).^(1/4)).*coeff,[0 0 4])
                        m=m+1
                    end

                    #H orbitals
                    if cmp(x2[l-1,1],"H")==0
                        #Hx5
                        arrayOforbitals[m]=orbitals("Gx5",expon,(((8388608*(expon.^13))./(893025*pi^3)).^(1/4)).*coeff,[5 0 0])
                        m=m+1
                        #Hx4y
                        arrayOforbitals[m]=orbitals("Gx4y",expon,(((8388608*(expon.^13))./(11025*pi^3)).^(1/4)).*coeff,[4 1 0])
                        m=m+1
                        #Hx4z
                        arrayOforbitals[m]=orbitals("Gx4y",expon,(((8388608*(expon.^13))./(11025*pi^3)).^(1/4)).*coeff,[4 0 1])
                        m=m+1
                        #Hx3y2
                        arrayOforbitals[m]=orbitals("Gx3y2",expon,(((8388608*(expon.^13))./(2025*pi^3)).^(1/4)).*coeff,[3 2 0])
                        m=m+1
                        #Hx2yz
                        arrayOforbitals[m]=orbitals("Gx4y",expon,(((8388608*(expon.^13))./(225*pi^3)).^(1/4)).*coeff,[3 1 1])
                        m=m+1
                        #Hx3z2
                        arrayOforbitals[m]=orbitals("Gx3z2",expon,(((8388608*(expon.^13))./(2025*pi^3)).^(1/4)).*coeff,[3 0 2])
                        m=m+1
                        #Hx2y3
                        arrayOforbitals[m]=orbitals("Gx2y3",expon,(((8388608*(expon.^13))./(2025*pi^3)).^(1/4)).*coeff,[2 3 0])
                        m=m+1
                        #Hx2y2z
                        arrayOforbitals[m]=orbitals("Gx2y2z",expon,(((8388608*(expon.^13))./(81*pi^3)).^(1/4)).*coeff,[2 2 1])
                        m=m+1
                        #Hx2yz2
                        arrayOforbitals[m]=orbitals("Gx2yz2",expon,(((8388608*(expon.^13))./(81*pi^3)).^(1/4)).*coeff,[2 1 2])
                        m=m+1
                        #Hx2z3
                        arrayOforbitals[m]=orbitals("Gx2z3",expon,(((8388608*(expon.^13))./(2025*pi^3)).^(1/4)).*coeff,[2 0 3])
                        m=m+1
                        #Hxy4
                        arrayOforbitals[m]=orbitals("Gxy4",expon,(((8388608*(expon.^13))./(11025*pi^3)).^(1/4)).*coeff,[1 4 0])
                        m=m+1
                        #Hxy3z
                        arrayOforbitals[m]=orbitals("Gxy3z",expon,(((8388608*(expon.^13))./(225*pi^3)).^(1/4)).*coeff,[1 3 1])
                        m=m+1
                        #Hxy2z2
                        arrayOforbitals[m]=orbitals("Gxy2z2",expon,(((8388608*(expon.^13))./(81*pi^3)).^(1/4)).*coeff,[1 2 2])
                        m=m+1
                        #Hxyz3
                        arrayOforbitals[m]=orbitals("Gxyz3",expon,(((8388608*(expon.^13))./(225*pi^3)).^(1/4)).*coeff,[1 1 3])
                        m=m+1
                        #Hxz4
                        arrayOforbitals[m]=orbitals("Gxz4",expon,(((8388608*(expon.^13))./(11025*pi^3)).^(1/4)).*coeff,[1 0 4])
                        m=m+1
                        #Hy5
                        arrayOforbitals[m]=orbitals("Gy5",expon,(((8388608*(expon.^13))./(893025*pi^3)).^(1/4)).*coeff,[0 5 0])
                        m=m+1
                        #Hy4z
                        arrayOforbitals[m]=orbitals("Gy4z",expon,(((8388608*(expon.^13))./(11025*pi^3)).^(1/4)).*coeff,[0 4 1])
                        m=m+1
                        #Hy3z2
                        arrayOforbitals[m]=orbitals("Gy3z2",expon,(((8388608*(expon.^13))./(2025*pi^3)).^(1/4)).*coeff,[0 3 2])
                        m=m+1
                        #Hy2z3
                        arrayOforbitals[m]=orbitals("Gy2z3",expon,(((8388608*(expon.^13))./(2025*pi^3)).^(1/4)).*coeff,[0 2 3])
                        m=m+1
                        #Hyz4
                        arrayOforbitals[m]=orbitals("Gyz4",expon,(((8388608*(expon.^13))./(11025*pi^3)).^(1/4)).*coeff,[0 1 4])
                        m=m+1
                        #Hz5
                        arrayOforbitals[m]=orbitals("Gz5",expon,(((8388608*(expon.^13))./(893025*pi^3)).^(1/4)).*coeff,[0 0 5])
                        m=m+1
                    end
                    l=l+x2[l-1,2]+1;
                end
                arrayOfAtoms[k]=atom(Geom,arrayOforbitals)
            end
        end
       return arrayOfAtoms,Mol
    end
