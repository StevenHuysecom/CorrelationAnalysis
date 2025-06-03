function LoadImages(file, chan)
    Folder = dir(file.path);
    %% Check if there was already a loading
    CheckForDir = [];
    for i = 3:size(Folder, 1)
        CheckForDir(end+1,1) = Folder(i).isdir;
    end

    if max(CheckForDir >= 1)
        disp('Data was already loaded and stored in folders')
    else
        %% make subfolder per sample       
        LifFiles = {};
        for i = 3:size(Folder, 1)
            if strcmp(Folder(i).name(end-2:end), 'lif')
                LifFiles{end+1,1} = Folder(i).name;
                subfolder = Folder(i).name(1:end-4);
                mkdir([file.path,filesep,subfolder])
    
                FileName = append([file.path,filesep,subfolder,file.ext]);
                bfI = BioformatsImage(FileName);
    
                 for k = 1:bfI.seriesCount
                    bfI.series = k;
                    ImageName = append('Position', num2str(k));
                    mkdir([file.path,filesep,subfolder,filesep,ImageName])
                    
                    PxSizes = bfI.pxSize;
                    MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep,'PxSizes.mat');
                    save(MatFileName, 'PxSizes')

                    for l = 1:bfI.sizeC
                        if l == 1
                            for j = 1:bfI.sizeZ
                                Membrane(:,:,j) = getPlane(bfI, j, l, 1);
                            end
                            MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep,'Membrane.mat');
                            save(MatFileName, 'Membrane')
                        elseif l == 2
                            for j = 1:bfI.sizeZ
                                Particles(:,:,j) = getPlane(bfI, j, l, 1);
                            end
                            MatFileName = append(file.path,filesep,subfolder,filesep,ImageName,filesep,'Particles.mat');
                            save(MatFileName, 'Particles')
                        else
                            assert("Unknown channel detected - not membrane neither particles")
                        end
                    end

                    clear Membrane Particles
                end
            end
        end
    end
end

