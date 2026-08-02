--!strict

local Floor = require(script.Parent:WaitForChild("Floor"))

local LocalGrid = {}

export type Cell = {
	ui: number, vi: number,   -- integer lattice indices in the part's local frame
	pos: Vector3,             -- exact surface position (world)
	normal: Vector3,
	slope: number,            -- degrees from world-up
	clearance: number,        -- studs of vertical headroom (capped)
	cover: Instance?,
}

export type DeadCell = {
	ui: number, vi: number,
	pos: Vector3,
	killer: Instance?,
}

export type Grid = {
	part: BasePart,
	fallback: boolean,        -- true => world-aligned (non-block part)
	origin: Vector3?,         -- face corner (world); block grids only
	u: Vector3?, v: Vector3?, -- in-plane unit axes (world); block grids only
	n: Vector3?,              -- surface normal (world); block grids only
	center: Vector3?,         -- centre of the walkable face (world)
	uExt: number?, vExt: number?, -- half-extents along u and v
	step: number,
	cells: {Cell},
	index: { [string]: Cell },-- "ui:vi" -> cell
	dead: {DeadCell},
	deadIndex: { [string]: DeadCell },
}

export type Config = {
	step: number?, maxSlope: number?, clearCap: number?, minClearance: number?,
}

local DEFAULT = {
	step = 1,           -- local cell size (studs)
	maxSlope = 65,      -- max walkable slope (deg); Cocosulx-tested
	clearCap = 20,      -- clearance raycast cap
	minClearance = 1.5, -- below this a cell isn't standable floor (crawl minimum)
}

local UP = Vector3.new(0, 1, 0)

local function merged(cfg): any
	local c = {}
	for k, v in pairs(DEFAULT) do c[k] = v end
	if cfg then for k, v in pairs(cfg) do if v ~= nil then c[k] = v end end end
	return c
end

local function isBlock(p: BasePart): boolean
	return p:IsA("Part") and p.Shape == Enum.PartType.Block
end

local function isClip(p: Instance): boolean
	return p.Name:find("ClipRamp") ~= nil
end

-- Group walkable surfels by the part beneath them.
local function groupByPart(surfels: {any}): { [BasePart]: {any} }
	local byPart: { [BasePart]: {any} } = {}
	for _, s in ipairs(surfels) do
		local b = byPart[s.part]
		if not b then b = {}; byPart[s.part] = b end
		b[#b + 1] = s
	end
	return byPart
end

-- Average the surfel normals to get the true walkable-face direction.
local function avgNormal(surfels: {any}): Vector3
	local s = Vector3.zero
	for _, sf in ipairs(surfels) do s += sf.normal end
	return (s.Magnitude > 1e-4) and s.Unit or Vector3.yAxis
end

local function topFace(part: BasePart, surfaceN: Vector3)
	local cf = part.CFrame
	local sz = part.Size
	local axes = {
		{ dir = cf.RightVector, ext = sz.X * 0.5 },
		{ dir = cf.UpVector,    ext = sz.Y * 0.5 },
		{ dir = cf.RightVector:Cross(cf.UpVector), ext = sz.Z * 0.5 }, -- local Z basis
	}
	local bi, best = 2, -math.huge
	for i, a in ipairs(axes) do
		local d = math.abs(a.dir:Dot(surfaceN))
		if d > best then best = d; bi = i end
	end
	local a = axes[bi]
	local n = (a.dir:Dot(surfaceN) >= 0) and a.dir or -a.dir
	local plane = {}
	for i, ax in ipairs(axes) do
		if i ~= bi then plane[#plane + 1] = ax end
	end
	return n, a.ext, plane[1], plane[2]
end

local function buildBlockGrid(part: BasePart, surfels: {any}, c: any, filterAll: RaycastParams, probe: BasePart, op: OverlapParams, rpTerrain: RaycastParams): Grid
	local n, nExt, ua, va = topFace(part, avgNormal(surfels))
	local u, uExt = ua.dir, ua.ext
	local v, vExt = va.dir, va.ext
	local surfaceCenter = part.Position + n * nExt
	local corner = surfaceCenter - u * uExt - v * vExt

	local rpPart = RaycastParams.new()
	rpPart.FilterType = Enum.RaycastFilterType.Include
	rpPart.FilterDescendantsInstances = { part }

	local grid: Grid = {
		part = part, fallback = false, origin = corner,
		u = u, v = v, n = n, step = c.step, cells = {}, index = {},
		dead = {}, deadIndex = {},
		center = surfaceCenter, uExt = uExt, vExt = vExt,
	}

	local function kill(iu: number, iv: number, pos: Vector3, killer: Instance?)
		local d: DeadCell = { ui = iu, vi = iv, pos = pos, killer = killer }
		grid.dead[#grid.dead + 1] = d
		grid.deadIndex[string.format("%d:%d", iu, iv)] = d
	end

	local step = c.step
	local nu = math.max(1, math.floor(2 * uExt / step + 1e-6))
	local nv = math.max(1, math.floor(2 * vExt / step + 1e-6))
	local castH = 2 -- studs above the surface to start the (downward-along-normal) ray

	for iu = 0, nu - 1 do
		for iv = 0, nv - 1 do
			local p = corner + u * ((iu + 0.5) * step) + v * ((iv + 0.5) * step)
			local res = workspace:Raycast(p + n * castH, -n * (castH + 0.5), rpPart)
			if not res then continue end
			local slope = math.deg(math.acos(math.clamp(res.Normal:Dot(UP), -1, 1)))
			if not ((slope <= c.maxSlope) or isClip(part)) then continue end
			probe.CFrame = CFrame.new(res.Position + UP * (0.1 + (c.minClearance - 0.1) * 0.5))
			local killer: Instance? = nil
			for _, hit in ipairs(workspace:GetPartsInPart(probe, op)) do
				if hit ~= part then killer = hit; break end
			end
			if killer then
				kill(iu, iv, res.Position, killer)
				continue
			end
			local upRes = workspace:Raycast(res.Position + Vector3.new(0, 0.15, 0), UP * c.clearCap, filterAll)
			local clearance = upRes and upRes.Distance or c.clearCap
			local cover: Instance? = upRes and upRes.Instance or nil
			local tUp = workspace:Raycast(res.Position + Vector3.new(0, 0.15, 0), UP * c.clearCap, rpTerrain)
			if tUp then
				if tUp.Distance < clearance then
					clearance = tUp.Distance
					cover = workspace.Terrain
				end
			elseif workspace:Raycast(res.Position + UP * c.clearCap, -UP * (c.clearCap - 0.25), rpTerrain) then
				clearance = 0
				cover = workspace.Terrain
			end
			if clearance < c.minClearance then
				kill(iu, iv, res.Position, cover)
				continue
			end
			local cell: Cell = {
				ui = iu, vi = iv, pos = res.Position, normal = res.Normal,
				slope = slope, clearance = clearance, cover = cover,
			}
			grid.cells[#grid.cells + 1] = cell
			grid.index[string.format("%d:%d", iu, iv)] = cell
		end
	end
	return grid
end

local function buildFallbackGrid(part: BasePart, surfels: {any}, c: any): Grid
	local grid: Grid = {
		part = part, fallback = true, step = c.step, cells = {}, index = {},
		dead = {}, deadIndex = {},
	}
	for _, s in ipairs(surfels) do
		if s.clearance < c.minClearance then continue end
		local iu = math.floor(s.pos.X / c.step)
		local iv = math.floor(s.pos.Z / c.step)
		local cell: Cell = {
			ui = iu, vi = iv, pos = s.pos, normal = s.normal,
			slope = s.slope, clearance = s.clearance, cover = s.cover,
		}
		grid.cells[#grid.cells + 1] = cell
		grid.index[string.format("%d:%d", iu, iv)] = cell
	end
	return grid
end

-- Build per-part local grids from an existing floor extraction.
function LocalGrid.fromFloor(floorData: any, parts: {BasePart}, cfg: Config?)
	local c = merged(cfg)
	local filterAll = RaycastParams.new()
	filterAll.FilterType = Enum.RaycastFilterType.Include
	filterAll.FilterDescendantsInstances = parts

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

	local byPart = groupByPart(floorData.surfels)
	local grids: { [BasePart]: Grid } = {}
	local nBlock, nFallback, nCells, nDead = 0, 0, 0, 0
	for part, sfs in pairs(byPart) do
		local g: Grid
		if isBlock(part) then
			g = buildBlockGrid(part, sfs, c, filterAll, probe, op, rpTerrain)
			nBlock += 1
		else
			g = buildFallbackGrid(part, sfs, c)
			nFallback += 1
		end
		grids[part] = g
		nCells += #g.cells
		nDead += #g.dead
	end
	probe:Destroy()

	return {
		grids = grids, config = c,
		stats = { parts = nBlock + nFallback, block = nBlock, fallback = nFallback, cells = nCells, dead = nDead },
	}
end

-- Convenience one-call bake: Floor.build + local grids.
-- Returns localData, floorData, tree, parts.
function LocalGrid.build(cfg: Config?)
	local floorData, tree, parts = Floor.build(cfg)
	local data = LocalGrid.fromFloor(floorData, parts, cfg)
	return data, floorData, tree, parts
end

function LocalGrid.visualize(data: any, parent: Instance?)
	local root = parent or workspace
	local dbg = root:FindFirstChild("NVGN_Debug")
	if not dbg then
		dbg = Instance.new("Folder"); dbg.Name = "NVGN_Debug"; dbg.Parent = root
	end
	local old = dbg:FindFirstChild("LocalGrid")
	if old then old:Destroy() end
	local folder = Instance.new("Folder"); folder.Name = "LocalGrid"; folder.Parent = dbg

	local step = data.config.step
	local i = 0
	for part, g in pairs(data.grids) do
		i += 1
		local hue = (i * 0.61803398875) % 1
		local sat = g.fallback and 0.3 or 0.9
		local pf = Instance.new("Folder"); pf.Name = part.Name; pf.Parent = folder
		for _, cell in ipairs(g.cells) do
			local dot = Instance.new("Part")
			dot.Anchored = true; dot.CanCollide = false; dot.CanQuery = false; dot.CanTouch = false
			local w, v
			if cell.clearance >= 4 then
				w, v = 0.9, 1
			elseif cell.clearance >= 3 then
				w, v = 0.7, 0.55
			else
				w, v = 0.55, 0.28
			end
			dot.Size = Vector3.new(w * step, 0.1, w * step)
			dot.Color = Color3.fromHSV(hue, sat, v)
			-- matte, so the neon Boundary edges pop over the grid layer
			dot.Material = Enum.Material.SmoothPlastic
			if not g.fallback and g.n then
				dot.CFrame = CFrame.fromMatrix(cell.pos, g.u, g.n)
			else
				dot.CFrame = CFrame.new(cell.pos)
			end
			dot.Name = string.format("c%.1f", cell.clearance)
			dot.Parent = pf
		end
	end
	return folder
end

return LocalGrid
