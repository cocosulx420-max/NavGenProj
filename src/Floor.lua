--!strict
-- NVGN.Floor — walkable surface extraction
--
-- Produces one surfel per 1-stud walkable cell. The SVO finds candidates
-- (solid voxels with empty space above); each candidate's top face is walked at
-- 1-stud resolution and a raycast onto the REAL part gives exact height + normal
-- (so ramps are smooth, not stair-stepped, and a collapsed node covering several
-- parts is sampled per-cell). Clearance = precise overlap probe (embedded-origin
-- detection) + upward raycast to the real ceiling; terrain via a separate ray pair.
-- Only walkable surfels are kept (steep faces feed the later boundary stage).

local SVO = require(script.Parent:WaitForChild("SVO"))

local Floor = {}

export type Surfel = {
	pos: Vector3,       -- exact surface position
	normal: Vector3,    -- surface normal
	slope: number,      -- degrees from world-up
	clearance: number,  -- studs of headroom above (capped; 0 = embedded/dead space)
	part: BasePart,     -- the part under this surfel
	-- A ClipRamp is the smooth surface authored over a staircase. Where one
	-- exists it IS the walkable surface and the steps beneath it are not; the
	-- boundary stage drops the risers wherever a clip surfel covers the cell.
	clip: boolean,
}

-- NOTE: horizontal "width" (distance to nearest wall) is intentionally NOT baked.
-- It is redundant with the navmesh boundary edges and is derived cheaply at
-- pathfinding time (portal-edge length + funnel radius offset). Clearance IS
-- baked because vertical headroom cannot be recovered from 2D boundaries.

export type Config = {
	leaf: number?, maxSlope: number?, agentHeight: number?,
	clearCap: number?, maxGroundFootprint: number?, minClearance: number?,
	root: Instance?, -- restrict the bake to this subtree (default: whole workspace)
}

local DEFAULT = {
	leaf = 1,                 -- SVO leaf size (studs)
	maxSlope = 65,            -- max walkable slope (deg); Cocosulx-tested
	agentHeight = 5,          -- reference stand height
	clearCap = 20,            -- clearance raycast cap
	maxGroundFootprint = 400, -- parts wider than this (baseplate) are excluded
	minClearance = 1.5,       -- headroom below this is dead space (crawl floor)
}

local UP = Vector3.new(0, 1, 0)

local function merged(cfg)
	local c = {}
	for k, v in pairs(DEFAULT) do c[k] = v end
	if cfg then for k, v in pairs(cfg) do if v ~= nil then c[k] = v end end end
	return c
end

local function isCharacter(p: Instance): boolean
	local a = p.Parent
	while a and a ~= workspace do
		if a:IsA("Model") and a:FindFirstChildOfClass("Humanoid") then return true end
		a = a.Parent
	end
	return false
end

local function cellKey(x: number, z: number): string
	return string.format("%d:%d", math.floor(x), math.floor(z))
end
Floor.cellKey = cellKey

-- Default world part filter: collidable, not terrain, not character,
-- not a huge flat ground slab (handled analytically elsewhere).
function Floor.gatherParts(cfg: Config?): {BasePart}
	local c = merged(cfg)
	local out = {}
	-- `root` scopes the bake to one model. A test map and the real map usually
	-- share a place, and baking the whole workspace to look at one of them wastes
	-- the only expensive stage in the pipeline.
	for _, d in ipairs((c.root or workspace):GetDescendants()) do
		if d:IsA("BasePart") and d.CanCollide and d.ClassName ~= "Terrain" and not isCharacter(d) then
			if math.max(d.Size.X, d.Size.Z) <= c.maxGroundFootprint then
				table.insert(out, d)
			end
		end
	end
	return out
end

-- Extract surfels from a prebuilt SVO over `parts`.
function Floor.extract(parts: {BasePart}, tree: any, cfg: Config?)
	local c = merged(cfg)
	local rp = RaycastParams.new()
	rp.FilterType = Enum.RaycastFilterType.Include
	rp.FilterDescendantsInstances = parts

	-- Embedded-origin probe: a thin invisible part spanning [0.1, minClearance]
	-- above each candidate, tested with precise GetPartsInPart (see clearance
	-- note below for why raycasts cannot do this job).
	local probe = Instance.new("Part")
	probe.Name = "NVGN_ClearProbe"
	probe.Size = Vector3.new(0.05, c.minClearance - 0.1, 0.05)
	probe.Anchored = true; probe.CanCollide = false; probe.CanQuery = false; probe.CanTouch = false
	probe.Transparency = 1
	probe.Parent = workspace
	local op = OverlapParams.new()
	op.FilterType = Enum.RaycastFilterType.Include
	op.FilterDescendantsInstances = parts
	local rpTerrain = RaycastParams.new()
	rpTerrain.FilterType = Enum.RaycastFilterType.Include
	rpTerrain.FilterDescendantsInstances = { workspace.Terrain }

	local surfels: {Surfel} = {}
	local index: { [string]: {Surfel} } = {}
	-- The same surface is now reachable from several leaves, so a column may be
	-- probed more than once. Key on cell and height so the index holds each
	-- surface once.
	local seen: {[string]: boolean} = {}

	-- The floor of the whole tree, not of this leaf. THIS IS THE POINT.
	--
	-- The down-ray used to be clipped to the leaf that spawned it (`edge + 2`),
	-- which silently assumed the surface under a column is inside the same leaf
	-- as whatever sits above it. At a wall it is not. The SVO is deliberately
	-- conservative, so a wall's voxels over-cover its footprint by up to a leaf;
	-- the ground leaf then skips that column on the empty-above guard (line
	-- below), and the wall's own leaf casts a ray far too short to reach the
	-- ground ten studs down. Neither leaf claims it, and the cell comes out with
	-- no surfel at all.
	--
	-- That is one missing cell hugging the base of every wall and the lip of
	-- every stair riser. Downstream they are not read as seams -- they are read
	-- as one-cell OBSTACLES, fall to the raw-ring fallback, and get drawn as a
	-- hexagon expanded by the agent radius. 100 of SmallMap's 102 "holes" were
	-- this, marching in a diagonal line down every edge in the map.
	--
	-- Casting to the bottom of the tree can only ADD surfels, never move one:
	-- the ray still stops at the first surface it meets. Clearance above them is
	-- measured against real geometry a few lines down, not against the
	-- conservative voxels, so a column that really is buried still reports
	-- clearance 0 and is filtered where it always was.
	local floorY = tree.center.Y - tree.half - 1

	tree:forEachSolidLeaf(function(ctr: Vector3, h: number)
		local edge = 2 * h
		local top = ctr.Y + h
		for i = 0, edge - 1 do
			for j = 0, edge - 1 do
				local cx = ctr.X - h + 0.5 + i
				local cz = ctr.Z - h + 0.5 + j
				if tree:isSolid(Vector3.new(cx, top + 0.5, cz)) then continue end
				local res = workspace:Raycast(Vector3.new(cx, top + 1, cz), Vector3.new(0, floorY - (top + 1), 0), rp)
				if not res then continue end
				local dk = cellKey(cx, cz) .. "|" .. math.floor(res.Position.Y * 4 + 0.5)
				if seen[dk] then continue end
				seen[dk] = true
				local n = res.Normal
				local slope = math.deg(math.acos(math.clamp(n:Dot(UP), -1, 1)))
				local isClip = res.Instance.Name:find("ClipRamp") ~= nil
				if not ((slope <= c.maxSlope) or isClip) then continue end
				-- Clearance. A raycast NEVER hits a part its origin is inside, so
				-- ray logic alone cannot detect an embedded origin (wall flush on
				-- the floor, buried overlap region, curved union). Precise overlap
				-- probe first: any foreign solid crossing [0.1, minClearance] above
				-- the surface means true headroom < minClearance -> clearance 0
				-- (surfel kept, truthful for the global index). Host excluded; a
				-- non-convex host's self-overhang is covered by the SVO empty-above
				-- guard at voxel scale. A clear probe guarantees the up-ray origin
				-- is outside every collider, making its distance exact.
				local clearance
				probe.CFrame = CFrame.new(res.Position + Vector3.new(0, 0.1 + (c.minClearance - 0.1) * 0.5, 0))
				local blocked = false
				for _, hit in ipairs(workspace:GetPartsInPart(probe, op)) do
					if hit ~= res.Instance then blocked = true; break end
				end
				if blocked then
					clearance = 0
				else
					local upRes = workspace:Raycast(res.Position + Vector3.new(0, 0.15, 0), Vector3.new(0, c.clearCap, 0), rp)
					clearance = upRes and upRes.Distance or c.clearCap
					-- Terrain is not in `parts` (never walkable) but still blocks
					-- headroom; overlap queries are parts-only, so use a terrain-only
					-- ray pair (down-ray catches embedded origin under a blob).
					local tUp = workspace:Raycast(res.Position + Vector3.new(0, 0.15, 0), Vector3.new(0, c.clearCap, 0), rpTerrain)
					if tUp then
						clearance = math.min(clearance, tUp.Distance)
					elseif workspace:Raycast(res.Position + Vector3.new(0, c.clearCap, 0), Vector3.new(0, -(c.clearCap - 0.25), 0), rpTerrain) then
						clearance = 0
					end
				end
				local surfel: Surfel = {
					pos = res.Position, normal = n, slope = slope,
					clearance = clearance, part = res.Instance, clip = isClip,
				}
				surfels[#surfels + 1] = surfel
				local key = cellKey(cx, cz)
				local bucket = index[key]
				if not bucket then bucket = {}; index[key] = bucket end
				bucket[#bucket + 1] = surfel
			end
		end
	end)
	probe:Destroy()

	return { surfels = surfels, index = index, config = c }
end

-- Convenience one-call bake: gather parts, build SVO, extract floor.
-- Returns floorData, tree, parts.
function Floor.build(cfg: Config?)
	local c = merged(cfg)
	local parts = Floor.gatherParts(c)
	local tree = SVO.fromParts(parts, c.leaf, 2)
	local data = Floor.extract(parts, tree, c)
	return data, tree, parts
end

return Floor
