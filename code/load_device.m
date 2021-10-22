function device = load_device
    load('../solution/deviceFile.mat', 'device');

    % Calculate r and xi for all measurement positions
    device.MEG.measurementPositionRadius = sqrt(sum(device.MEG.measurementPosition.^2,2));
    device.MEG.measurementPositionNormalized = device.MEG.measurementPosition./repmat(device.MEG.measurementPositionRadius,1,3);
end