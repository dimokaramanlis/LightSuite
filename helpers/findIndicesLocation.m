function locations = findIndicesLocation(indices, array)
  % Finds the locations of a group of indices in an array.
  %
  % Args:
  %   indices: A vector of values to search for in the array.
  %   array:   The array to search within.
  %
  % Returns:
  %   locations: A vector containing the locations (indices) of the 'indices'
  %              values within the 'array'.  Returns an empty array if any
  %              of the 'indices' are not found.  Handles duplicates correctly.
  %              Returns the *first* occurrence if a value in `indices` is 
  %              duplicated in `array`.
  %

  locations = zeros(1, length(indices)); % Pre-allocate for efficiency

  for i = 1:length(indices)
    index = find(array == indices(i), 1, 'first'); % Find the *first* occurrence
    if isempty(index)
      locations = []; % Return empty array if any index is not found
      return;
    end
    locations(i) = index;
  end
end
